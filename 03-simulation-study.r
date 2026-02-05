########################################
# Salvage therapy simulation study
# Lukas Owens 2025/12/1
#
# Goal: Create simulated longitudinal biomarker and event data where treatment
# initiation is driven by biomarker. Fit models across each simulated dataset
# and assess validity of estimation process.
########################################
library(tidyverse)
library(JMbayes2)
library(nlme)

# Output directory
out_dir <- 'output'

# Run this script N times with command line arguments 1 through N=500 to create
# N simulated datasets and models
index <- as.integer(commandArgs(trailingOnly = TRUE)[1])

########################################
# Simulation parameters
########################################
# Number of subjects
N <- 2000
# Number of measurements per subject
K <- 25
# Maximum follow-up time
t_max <- 8

########################################
# Model structure formulas
########################################
formula_beta <- ~ 1 +
    w1 +
    w2 +
    time +
    time:w1 +
    time:w2

formula_b <- ~ 1 + time

formula_gamma_mets <- ~ -1 +
    w1 +
    w2 +
    salvage_indicator +
    salvage_indicator:psa_value_st

formula_alpha_mets <- ~ -1 +
    value_logpsa:as.numeric(salvage_indicator==0)

########################################
# Model parameters
########################################
beta <- c('Intercept'=0.125,
          'w1'=0.04,
          'w2'=-0.02,
          'time'=0.150,
          'w1:time'=0.30,
          'w2:time'=-0.030)
D <- matrix(c(0.01, 0.03,
              0.03, 0.15),
            ncol=2)
sigma <- 0.05

# Weibull parameters for events
shapes <- c('mets'=3, 'death'=1)
scales <- c('mets'=10, 'death'=20)

# Hazard model parameters
gamma <- c('w1'=0.5,
           'w2'=-0.1,
           'salvage_indicator'=-0.5,
           'salvage_indicator:psa_value_st'=0.8)
alpha <- c('value_logpsa:!salvage_indicator'=1.0)

########################################
# Helper functions
########################################
# Probability of starting ST at a given observed PSA level
st_prob <- function(psa) {
    case_when(psa<0.01 ~ 0.01,
              psa < 1.0 ~ 0.15,
              psa < 2.0 ~ 0.20,
              TRUE ~ 0.4)
}

# Simulate salvage therapy time for an individual. A specified fraction of
# patients will never start therapy. For the remaining patients, treatment will
# be started with probability based on the `st_prob` function
simulate_salvage_time <- function(pset, pr_never=0.20) {
    pr <- st_prob(pset$psa)
    A <- runif(pr) < pr
    never <- runif(1) < pr_never
    if(all(!A) | never)
        return(100)
    else
        return(pset$time[which.max(A)])
}

# Return a function for the mean PSA growth based on fixed and random effects
eta_factory <- function(w1, w2, b) {
    function(time) {
        dset <- tibble(time=time, w1=w1, w2=w2, b=matrix(b, ncol=2))
        fixed_eff <- as.vector(model.matrix(formula_beta, dset) %*% beta)
        random_eff <- rowSums(model.matrix(formula_b, dset) %*% matrix(b, ncol=1))
        as.vector(fixed_eff + random_eff)
    }
}
# Usage: 
# f <- eta_factory(0, 2, c(0.1, 0.3))
# f(c(1,2,3))

# Create observed PSAs
simulate_long_data <- function(w1, w2, b) {
    pset <- tibble(time=cumsum(rexp(K, 3)))
    eta <- eta_factory(w1, w2, b)
    pset <- pset |> mutate(logpsa=eta(time))
    pset <- pset |> mutate(logpsa=logpsa+rnorm(nrow(pset), 0, sigma))
    pset <- pset |> mutate(psa=exp(logpsa)-1)
    pset
}

# Function to invert to get simulated time of metastasis
invert_this <- function(time, dset, i, event='mets') {
    shape <- shapes[event]
    scale <- scales[event]

    tset <- dset[i,]
    eta <- eta_factory(tset$w1, tset$w2, tset$b)

    h <- function(s) {
        tset <- bind_cols(tset, tibble(time=s))
        tset <- tset |> mutate(salvage_indicator=as.integer(time >= salvage_time))
        tset <- tset |> mutate(value_logpsa=eta(time))
        basehaz <- (shape/scale)*(s/scale)^(shape-1)
        loghr <- model.matrix(formula_gamma_mets, tset) %*% gamma + model.matrix(formula_alpha_mets, tset) %*% alpha
        basehaz*exp(loghr)
    }
    integrate(h, 0, time)$value + log(tset$u)
}

########################################
# Create all simulated patients
########################################
make_all_data <- function(seed) {
    # Simulate baseline covariates and PSA measurements
    set.seed(seed)
    dset <- tibble(id=1:N)
    dset <- dset |> mutate(w1=sample(0:1, size=N, replace=TRUE, prob=(2:1)/3))
    dset <- dset |> mutate(w2=rexp(N, 0.2))
    dset <- dset |> mutate(b=map(id, \(x) MASS::mvrnorm(1, mu=c(0,0), Sigma=D)))
    dset <- dset |> mutate(b=MASS::mvrnorm(N, mu=c(0,0), Sigma=D))
    dset <- dset |> rowwise()
    dset <- dset |> mutate(long_data=list(simulate_long_data(w1, w2, b)))

    # Simulate salvage therapy times
    dset$salvage_time <- map_dbl(dset$long_data, simulate_salvage_time)

    # Calculate eta and eta' at salvage
    dset <- dset |> rowwise()
    dset <- dset |> mutate(psa_value_st=map_dbl(salvage_time, eta_factory(w1, w2, b)))
    dset <- dset |> ungroup()

    # Simulate event times
    dset <- dset |> mutate(u=runif(nrow(dset)))
    f <- function(i) {
        x <- try(uniroot(invert_this,
                             interval=c(1e-5, 50),
                             dset=dset,
                             i=i)$root, TRUE)
        if(inherits(x, 'try-error')) return(50) else return(x)
    }

    times_mets <- map_dbl(1:N, f)
    err_rate <- sum(times_mets==50)
    cat('Number mets greater than 50: ', err_rate, '\n')
    times_death <- rweibull(N, shape=shapes['death'], scale=scales['death'])
    times_censor <- rep(t_max, N)
    dset$event <- case_when(times_censor < pmin(times_death, times_mets) ~ 0,
                            times_mets < pmin(times_death, times_censor) ~ 1,
                            times_death < pmin(times_mets, times_censor) ~ 2)
    dset <- dset |> mutate(time=case_when(event==0 ~ times_censor,
                                          event==1 ~ times_mets,
                                          event==2 ~ times_death))
    save(dset, file=str_glue('{out_dir}/data_all_{seed}.RData'))
}

########################################
# Fit model to simulated data
########################################
dummy <- function (f, lvl) as.numeric(f == lvl)
fit_models <- function(seed) {
    # Load data for selected seed
    load(str_glue('{out_dir}/data_all_{seed}.RData'))

    # Prep longitudinal data
    data_long <- dset |> select(id, w1, w2, long_data, salvage_time)
    data_long <- data_long |> unnest(long_data)

    # Estimate initial longitudinal model
    mod_long <- lme(data = data_long,
                    control = nlme::lmeControl(opt = 'optim'),
                    method = 'ML',
                    fixed = logpsa ~ 1 +
                        w1 +
                        w2 +
                        time +
                        w1:time +
                        w2:time,
                    random = list(id = pdDiag(form = ~ 1 + time)))

    # Prep survival data
    # Attach mean PSA at salvage
    data_surv <- dset |> select(id, w1, w2, salvage_time, event, time)
    beta_hat <- mod_long$coefficients$fixed
    m0 <- model.matrix(~w1 + w2 + salvage_time + w1:salvage_time + w2:salvage_time, data_surv)
    data_surv$psa_value_st  <- as.vector(m0 %*% beta_hat +
                                         mod_long$coefficients$random$id[,1] +
                                         mod_long$coefficients$random$id[,2] * data_surv$salvage_time)

    # Format data for survival analysis
    data_surv <- data_surv |> mutate(time_1=ifelse(salvage_time<time, salvage_time, time),
                                     time_2=ifelse(salvage_time<time, time, NA),
                                     event_1=ifelse(salvage_time<time, 0, event),
                                     event_2=ifelse(salvage_time<time, event, NA),
                                     salvage_1=0,
                                     salvage_2=1)

    data_surv <- data_surv |> select(-time, -event)
    data_surv <- data_surv |> pivot_longer(cols=starts_with(c('time', 'event', 'salvage')),
                                           names_to=c('.value', 'index'),
                                           names_sep='_')
    data_surv <- data_surv |> filter(!is.na(time))
    data_surv <- data_surv |> group_by(id)
    data_surv <- data_surv |> mutate(start=lag(time, default=0),
                                     stop=time)
    data_surv <- crossing(data_surv, tibble(cr=c('mets', 'death')))
    data_surv <- data_surv |> mutate(event=case_when(event==0~0,
                                                     event==1 & cr=='mets' ~ 1,
                                                     event==2 & cr=='death' ~ 1,
                                                     TRUE ~ 0))

    # Estimate initial survival model
    dummy <- function (f, lvl) as.numeric(f == lvl)
    mod_surv <- coxph(Surv(start, stop, event) ~
                      strata(cr) +
                      dummy(cr, 'mets'):w1 +
                      dummy(cr, 'mets'):w2 +
                      dummy(cr, 'mets'):salvage +
                      dummy(cr, 'mets'):salvage:psa_value_st +
                      dummy(cr, 'death'):w1,
                  data=data_surv)

    # Estimate joint model
    mod_jm <- jm(mod_surv,
                 mod_long,
                 data_Surv = data_surv,
                 data_Long = data_long,
                 time_var = 'time',
                 n_iter = 5000,
                 n_burnin = 100,
                 functional_forms = ~ dummy(cr, 'mets'):dummy(salvage, 0):value(logpsa) +
                     dummy(cr, 'death'):value(logpsa))

    # Save outputs
    save(data_long, file=str_glue('{out_dir}/data_long_{seed}.RData'))
    save(data_surv, file=str_glue('{out_dir}/data_surv_{seed}.RData'))
    save(mod_long, file=str_glue('{out_dir}/mod_long_{seed}.RData'))
    save(mod_surv, file=str_glue('{out_dir}/mod_surv_{seed}.RData'))
    save(mod_jm, file=str_glue('{out_dir}/mod_jm_{seed}.RData'))
}

########################################
# Run analysis
########################################
make_all_data(index)
fit_models(index)

