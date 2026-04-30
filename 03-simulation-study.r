########################################
# Simulation study for predictimands in
# salvage therapy joint models
# Lukas Owens 2025/12/1
########################################
library(tidyverse)
library(JMbayes2)
library(nlme)
library(future)
library(furrr)

plan(multisession)
options(future.rng.onMisuse='ignore')

out_dir <- 'output'
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

########################################
# Simulation parameters
########################################
N_sim <- 500          # simulation replications
N <- 2000             # subjects per simulated dataset
K <- 25               # PSA measurements per subject
t_max <- 8            # maximum follow-up (years)

N_patients <- 50      # patients for predictimand evaluation
n_long_samples <- 20  # simulated PSA histories per patient
n_inner_samples <- 30 # posterior draws per estimated predictimand
time_now <- 1
time_end <- time_now + 5

########################################
# True model parameters
########################################
theta_true <- list(
    beta=c('Intercept'=0.125,
           'w1'=0.04,
           'w2'=-0.02,
           'time'=0.150,
           'w1:time'=0.30,
           'w2:time'=-0.030),
    D=matrix(c(0.01, 0.03,
               0.03, 0.15),
             ncol=2),
    sigma=0.05,
    gamma=c('w1'=0.5,
            'w2'=-0.1,
            'salvage_indicator'=-0.5,
            'salvage_indicator:psa_value_st'=0.8),
    alpha=c('value_logpsa:!salvage_indicator'=1.0)
)
shapes <- c('mets'=3, 'death'=1)
scales <- c('mets'=10, 'death'=20)

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
# Data generation functions
########################################
st_prob <- function(psa) {
    case_when(psa < 0.01 ~ 0.01,
              psa < 1.0  ~ 0.15,
              psa < 2.0  ~ 0.20,
              TRUE       ~ 0.4)
}

simulate_salvage_time <- function(pset, pr_never=0.20) {
    pr <- st_prob(pset$psa)
    A <- runif(length(pr)) < pr
    never <- runif(1) < pr_never
    if(all(!A) | never)
        return(100)
    else
        return(pset$time[which.max(A)])
}

eta_factory <- function(w1, w2, b, beta=theta_true$beta) {
    function(time) {
        dset <- tibble(time=time, w1=w1, w2=w2)
        fixed_eff <- as.vector(model.matrix(formula_beta, dset) %*% beta)
        random_eff <- rowSums(model.matrix(formula_b, dset) %*% matrix(b, ncol=1))
        as.vector(fixed_eff + random_eff)
    }
}

# Simulate real-world PSA data at irregular intervals
simulate_long_data <- function(w1, w2, b) {
    pset <- tibble(time=cumsum(rexp(K, 3)))
    eta <- eta_factory(w1, w2, b)
    pset <- pset |> mutate(logpsa=eta(time))
    pset <- pset |> mutate(logpsa=logpsa + rnorm(nrow(pset), 0, theta_true$sigma))
    pset <- pset |> mutate(psa=exp(logpsa) - 1)
    pset
}

# Simulate PSAs at regular intervals when following specific strategy
simulate_psa_grid <- function(w1, w2, b, time_stop) {
    pset <- tibble(time=seq(0, time_stop, by=1/3))
    eta <- eta_factory(w1, w2, b)
    pset <- pset |> mutate(logpsa=eta(time))
    pset <- pset |> mutate(logpsa=logpsa + rnorm(nrow(pset), 0, theta_true$sigma))
    pset <- pset |> mutate(psa=exp(logpsa) - 1)
    pset
}


# Simulate time of metastasis given baseline covariates and random effects
invert_this <- function(time, u, w1, w2, b, salvage_time=100) {
    shape <- shapes['mets']
    scale <- scales['mets']
    eta <- eta_factory(w1, w2, b)
    tset <- tibble(w1=w1, w2=w2, salvage_time=salvage_time, psa_value_st=eta(salvage_time))
    h <- function(s) {
        tset <- bind_cols(tset, tibble(time=s))
        tset <- tset |> mutate(salvage_indicator=as.integer(time >= salvage_time))
        tset <- tset |> mutate(value_logpsa=eta(time))
        basehaz <- (shape/scale)*(s/scale)^(shape-1)
        loghr <- model.matrix(formula_gamma_mets, tset) %*% theta_true$gamma +
                 model.matrix(formula_alpha_mets, tset) %*% theta_true$alpha
        basehaz*exp(loghr)
    }
    integrate(h, 0, time)$value + log(u)
}
simulate_time_mets <- function(w1, w2, b, salvage_time=100) {
    u <- runif(1)
    x <- try(uniroot(invert_this,
                     interval=c(1e-5, 50),
                     u=u, w1=w1, w2=w2, b=b,
                     salvage_time=salvage_time)$root, TRUE)
    if(inherits(x, 'try-error')) return(50) else return(x)
}

dummy <- function(f, lvl) as.numeric(f == lvl)

########################################
# Generate simulated dataset
########################################
make_all_data <- function(seed) {
    set.seed(seed)
    dset <- tibble(id=1:N)
    dset <- dset |> mutate(w1=sample(0:1, size=N, replace=TRUE, prob=(2:1)/3))
    dset <- dset |> mutate(w2=rexp(N, 0.2))
    dset <- dset |> mutate(b=MASS::mvrnorm(N, mu=c(0,0), Sigma=theta_true$D))
    dset <- dset |> rowwise()
    dset <- dset |> mutate(long_data=list(simulate_long_data(w1, w2, b)))
    dset$salvage_time <- map_dbl(dset$long_data, simulate_salvage_time)
    dset <- dset |> ungroup()
    times_mets <- map_dbl(1:N, \(i) simulate_time_mets(dset$w1[i], dset$w2[i],
                                                        dset$b[i,], dset$salvage_time[i]))
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
# Fit joint model to simulated data
########################################
fit_models <- function(seed) {
    load(str_glue('{out_dir}/data_all_{seed}.RData'))

    data_long <- dset |> select(id, w1, w2, long_data, salvage_time, event_time=time)
    data_long <- data_long |> unnest(long_data)
    data_long <- data_long |> filter(time <= salvage_time, time <= event_time)

    mod_long <- lme(data=data_long,
                    control=nlme::lmeControl(opt='optim'),
                    method='ML',
                    fixed=logpsa ~ 1 + w1 + w2 + time + w1:time + w2:time,
                    random=list(id=pdDiag(form=~ 1 + time)))

    data_surv <- dset |> select(id, w1, w2, salvage_time, event, time)
    data_surv <- data_surv |> filter(id %in% unique(data_long$id))
    beta_hat <- mod_long$coefficients$fixed
    m0 <- model.matrix(~w1 + w2 + salvage_time + w1:salvage_time + w2:salvage_time, data_surv)
    data_surv$psa_value_st <- as.vector(m0 %*% beta_hat +
                                        mod_long$coefficients$random$id[,1] +
                                        mod_long$coefficients$random$id[,2] * data_surv$salvage_time)

    data_surv <- data_surv |> mutate(time_1=ifelse(salvage_time < time, salvage_time, time),
                                     time_2=ifelse(salvage_time < time, time, NA),
                                     event_1=ifelse(salvage_time < time, 0, event),
                                     event_2=ifelse(salvage_time < time, event, NA),
                                     salvage_1=0,
                                     salvage_2=1)
    data_surv <- data_surv |> select(-time, -event)
    data_surv <- data_surv |> pivot_longer(cols=starts_with(c('time', 'event', 'salvage')),
                                           names_to=c('.value', 'index'),
                                           names_sep='_')
    data_surv <- data_surv |> filter(!is.na(time))
    data_surv <- data_surv |> group_by(id)
    data_surv <- data_surv |> mutate(start=lag(time, default=0), stop=time)
    data_surv <- crossing(data_surv, tibble(cr=c('mets', 'death')))
    data_surv <- data_surv |> mutate(event=case_when(event==0 ~ 0,
                                                     event==1 & cr=='mets' ~ 1,
                                                     event==2 & cr=='death' ~ 1,
                                                     TRUE ~ 0))

    mod_surv <- coxph(Surv(start, stop, event) ~
                          strata(cr) +
                          dummy(cr, 'mets'):w1 +
                          dummy(cr, 'mets'):w2 +
                          dummy(cr, 'mets'):salvage +
                          dummy(cr, 'mets'):salvage:psa_value_st +
                          dummy(cr, 'death'):w1,
                      data=data_surv)

    mod_jm <- jm(mod_surv,
                 mod_long,
                 data_Surv=data_surv,
                 data_Long=data_long,
                 time_var='time',
                 n_iter=5000,
                 n_burnin=100,
                 functional_forms=~ dummy(cr, 'mets'):dummy(salvage, 0):value(logpsa) +
                     dummy(cr, 'death'):value(logpsa),
                 save_random_effects=FALSE)

    save(data_long, file=str_glue('{out_dir}/data_long_{seed}.RData'))
    save(data_surv, file=str_glue('{out_dir}/data_surv_{seed}.RData'))
    save(mod_long, file=str_glue('{out_dir}/mod_long_{seed}.RData'))
    save(mod_surv, file=str_glue('{out_dir}/mod_surv_{seed}.RData'))
    save(mod_jm, file=str_glue('{out_dir}/mod_jm_{seed}.RData'))
}

########################################
# Predictimand utilities
########################################
# Integrate a function with a new discontinuity (faster than naive approach)
integrate_discont <- function(f, lower, upper, discont=NULL) {
    if (!is.null(discont) & lower < discont & discont < upper) {
        return(integrate(f, lower, discont)$value + integrate(f, discont, upper)$value)
    } else {
        return(integrate(f, lower, upper)$value)
    }
}

# Log-likelihood of observed data for sampling procedure
loglik <- function(b, w1, w2, theta, long_data, time_now) {
    eta <- eta_factory(w1, w2, b)
    l_long <- -1/theta$sigma^2 * sum((long_data$logpsa - eta(long_data$time))^2)
    tset <- tibble(w1=w1, w2=w2, salvage_indicator=0, psa_value_st=0)
    h <- function(s) {
        tset <- bind_cols(tset, tibble(time=s))
        tset <- tset |> mutate(value_logpsa=eta(time))
        basehaz_mets <- (shapes['mets']/scales['mets'])*(s/scales['mets'])^(shapes['mets']-1)
        basehaz_death <- (shapes['death']/scales['death'])*(s/scales['death'])^(shapes['death']-1)
        loghr <- model.matrix(formula_gamma_mets, tset) %*% theta$gamma +
                 model.matrix(formula_alpha_mets, tset) %*% theta$alpha
        basehaz_mets*exp(loghr) + basehaz_death
    }
    l_surv <- -integrate(h, lower=0, upper=time_now)$value
    l_prior <- -1/2 * as.vector(b %*% solve(theta$D) %*% t(b))
    l_long + l_surv + l_prior
}

# Sample random effects using Metropolis-Hastings
sample_b <- function(log_target,
                     b_init=matrix(c(0,0), nrow=1),
                     n_iter=1200,
                     prop_scale=c(0.001, 0.001),
                     burnin=floor(0.20*n_iter),
                     thin=1) {
    d <- length(b_init)
    prop_cov <- diag(prop_scale, d)
    b_current <- b_init
    loglik_current <- log_target(b_current)
    draws <- matrix(NA_real_, nrow=n_iter, ncol=d)
    accept <- 0L
    for (iter in seq_len(n_iter)) {
        b_next <- matrix(MASS::mvrnorm(1, mu=b_current, Sigma=prop_cov), nrow=1)
        loglik_next <- log_target(b_next)
        log_alpha <- loglik_next - loglik_current
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            b_current <- b_next
            loglik_current <- loglik_next
            accept <- accept + 1L
        }
        draws[iter, ] <- b_current
    }
    keep <- seq.int(from=burnin+1, to=n_iter, by=thin)
    list(draws=draws[keep, , drop=FALSE],
         accept_rate=accept/n_iter)
}

# Generate treatment time according to a specified strategy
treatment_time <- function(strategy, w1, w2, b, time_now, time_end, threshold=2) {
    if(strategy == 'now') return(time_now)
    if(strategy == 'never') return(100)
    if(strategy == 'dynamic') {
        psa_times <- seq(time_now, time_end, by=1/3)
        eta <- eta_factory(w1, w2, b)
        logpsas <- eta(psa_times) + rnorm(length(psa_times), 0, theta_true$sigma)
        psas <- exp(logpsas) - 1
        if(all(psas < threshold))
            return(100)
        else
            return(psa_times[which.max(psas >= threshold)])
    }
}

########################################
# Compute predictimands for one fitted model
########################################
run_predictimands <- function(seed) {
    set.seed(seed)
    load(str_glue('{out_dir}/mod_jm_{seed}.RData'))

    bs_knots <- mod_jm$control$knots[[1]]
    bs_ord <- mod_jm$control$Bsplines_degree[1] + 1
    n_bs_per_risk <- ncol(mod_jm$mcmc$bs_gammas[[1]]) / 2
    theta_sample <- list(
        gammas=mod_jm$mcmc$gammas[[1]],
        alphas=mod_jm$mcmc$alphas[[1]],
        beta=mod_jm$mcmc$betas1[[1]],
        bs_gammas=mod_jm$mcmc$bs_gammas[[1]],
        W_std_gammas=mod_jm$mcmc$W_std_gammas[[1]],
        Wlong_std_alphas=mod_jm$mcmc$Wlong_std_alphas[[1]],
        length=mod_jm$control$n_iter - mod_jm$control$n_burnin
    )

    # Baseline hazard using spline from fitted model.
    bh <- function(t, theta, risk='mets') {
        idx <- if (risk == 'mets') (n_bs_per_risk + 1):(2 * n_bs_per_risk) else 1:n_bs_per_risk
        W0 <- splineDesign(knots=bs_knots, x=t, ord=bs_ord, outer.ok=TRUE)
        log_bh <- c(W0 %*% theta$bs_gammas[idx]) - theta$W_std_gammas - theta$Wlong_std_alphas
        exp(log_bh)
    }

    sample_theta <- function(s) {
        set.seed(s)
        idx <- sample(1:theta_sample$length, size=1)
        gamma_names <- c('dummy(cr, "mets"):w1',
                         'dummy(cr, "mets"):w2',
                         'dummy(cr, "mets"):salvage',
                         'dummy(cr, "mets"):salvage:psa_value_st')
        beta_names <- c('(Intercept)', 'w1', 'w2', 'time', 'w1:time', 'w2:time')
        list(gamma=theta_sample$gammas[idx, gamma_names],
             alpha=theta_sample$alphas[idx, 'value(logpsa):dummy(cr, "mets"):dummy(salvage, 0)'],
             beta=theta_sample$beta[idx, beta_names],
             bs_gammas=theta_sample$bs_gammas[idx, ],
             W_std_gammas=theta_sample$W_std_gammas,
             Wlong_std_alphas=theta_sample$Wlong_std_alphas)
    }

    pr_func <- function(time, w1, w2, b, theta, salvage_time) {
        use_spline <- !is.null(theta$bs_gammas)
        eta <- eta_factory(w1, w2, b, beta=theta$beta)
        dset <- tibble(w1=w1, w2=w2, salvage_time=salvage_time, psa_value_st=eta(salvage_time))
        hE <- function(s) {
            dset <- bind_cols(dset, tibble(time=s))
            dset <- dset |> mutate(salvage_indicator=as.integer(time >= salvage_time))
            dset <- dset |> mutate(value_logpsa=eta(time))
            basehaz <- if (use_spline) {
                bh(s, theta, risk='mets')
            } else {
                (shapes['mets']/scales['mets'])*(s/scales['mets'])^(shapes['mets']-1)
            }
            loghr <- model.matrix(formula_gamma_mets, dset) %*% theta$gamma +
                     model.matrix(formula_alpha_mets, dset) %*% theta$alpha
            basehaz*exp(loghr)
        }
        hD <- function(s) {
            if (use_spline) {
                bh(s, theta, risk='death')
            } else {
                (shapes['death']/scales['death'])*(s/scales['death'])^(shapes['death']-1)
            }
        }
        surv <- function(s) {
            exp(-integrate_discont(\(u) hE(u) + hD(u), lower=0, upper=s, discont=salvage_time))
        }
        ts <- seq(time_now, time, by=0.1)
        ys <- map_dbl(ts, surv)
        surv_interp <- approxfun(x=ts, y=ys)
        f <- function(s) hE(s) * surv_interp(s) / surv_interp(time_now)
        integrate_discont(f, time_now, time, discont=salvage_time)
    }

    predictimand <- function(w1, w2, b, strategy, theta) {
        salvage_time <- treatment_time(strategy, w1, w2, b, time_now, time_end)
        pr_func(time_end, w1, w2, b, theta, salvage_time)
    }

    predictimand_estimated <- function(w1, w2, b_samples, strategy) {
        map_dbl(1:n_inner_samples, \(s) {
            theta <- sample_theta(s)
            set.seed(s)
            b <- b_samples[sample(1:nrow(b_samples), 1), ]
            predictimand(w1, w2, b, strategy, theta)
        })
    }

    # Set up simulated patient cohort
    dset <- tibble(id=1:N_patients)
    dset <- dset |> mutate(w1=sample(0:1, size=N_patients, replace=TRUE, prob=(2:1)/3))
    dset <- dset |> mutate(w2=rexp(N_patients, 0.2))
    dset <- dset |> mutate(b=MASS::mvrnorm(N_patients, mu=c(0,0), Sigma=theta_true$D))
    dset <- dset |> rowwise()
    dset <- dset |> mutate(time_mets=simulate_time_mets(w1, w2, b))
    dset <- dset |> filter(time_mets > time_now)
    dset <- dset |> crossing(long_data_idx=1:n_long_samples)
    dset <- dset |> rowwise()
    dset <- dset |> mutate(long_data=list(simulate_psa_grid(w1, w2, b, time_now)))

    # Sample posterior random effects given observed PSA history
    dset$b_samples <- pmap(list(dset$w1, dset$w2, dset$long_data), \(w1, w2, long_data) {
        target <- \(x) loglik(x, w1, w2, theta_true, long_data, time_now)
        sample_b(target)$draws
    })

    # Compute true and estimated predictimands under each strategy
    dset <- dset |> crossing(strategy=c('now', 'never', 'dynamic'))
    b_list <- map(1:nrow(dset), \(x) dset$b[x,])
    dset$pred_true <- pmap_dbl(list(dset$w1, dset$w2, b_list, dset$strategy),
                               \(w1, w2, b, strategy) predictimand(w1, w2, b, strategy, theta_true))
    dset$pred_est <- pmap(list(dset$w1, dset$w2, dset$b_samples, dset$strategy),
                          \(w1, w2, b_samples, strategy) predictimand_estimated(w1, w2, b_samples, strategy))

    save(dset, file=str_glue('{out_dir}/predictimands-out-{seed}.RData'))
}

########################################
# Step 1: Generate data and fit models
########################################
future_walk(1:N_sim, \(seed) {
    make_all_data(seed)
    fit_models(seed)
})

########################################
# Step 2: Compute predictimands
########################################
future_walk(1:N_sim, run_predictimands)
