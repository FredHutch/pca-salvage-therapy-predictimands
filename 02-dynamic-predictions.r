###################################
# 02-dynamic-predictions.R
#
# Dynamic predictions of patient outcomes under alternative salvage therapy
# treatment strategies
#
# PSA: Prostate-specific antigen
# ST: Salvage therapy
# RP: Radical prostatectomy
# BCR: Biochemical recurrence
###################################

library(tidyverse)
library(JMbayes2)
library(furrr)
set.seed(1234)

########################################
# Load coefficients and specify model structure
########################################
# Which functional form to use (1: value only, 2: slope only, 3: both)
ff <- 3

# Load fitted model
load(str_glue('02-fit-model-1/output/model-jm-ff{ff}.RData'))

# Helper functions for factor variables
cci_fct <- \(x) factor(x, levels = c('0', '1', '2+'))
tx_grade_group_fct <- \(x) factor(x, levels = c('<=6', '3+4', '4+3', '8', '9-10'))
tx_tstage_fct <- \(x) factor(x, levels = c('T0-T2', 'T3a', 'T3b', 'T4'))
sms_lni_fct <- \(x) factor(x, levels = c('Negative', 'Positive'))
dummy <- function (f, lvl) as.numeric(f == lvl)

# Model formulas
formula_beta <- ~ 1 +
    tx_grade_group + 
    lni + 
    I(dx_logpsa - 2) + 
    tt_bcr + 
    psa_time + 
    I(psa_time^2) + 
    psa_time:tx_grade_group + 
    psa_time:lni + 
    psa_time:tt_bcr

formula_b <- ~ 1 + 
    psa_time + 
    I(psa_time^2)

# Formulas for constant, linear and quadratic time terms
formula_beta_0 <- ~ 1 + tx_grade_group + lni + I(dx_logpsa-2) + tt_bcr
formula_beta_1 <- ~ 1 + tx_grade_group + lni + tt_bcr
formula_beta_2 <- ~ 1

if (ff==1) {
    formula_gamma_mets <- ~ 1 +
        I(bcr_age - 50) +
        I(dx_logpsa - 2) + 
        tx_grade_group + 
        tx_tstage + 
        sms + 
        lni + 
        tt_bcr + 
        st_ht + 
        st_rt + 
        value_logpsa_st:salvage_indicator +
        st_ht:st_rt
    formula_alpha_mets <- ~ -1 +
        value_logpsa:as.numeric(salvage_indicator==0)
    
} else if (ff==2) {
    formula_gamma_mets <- ~ 1 +
        I(bcr_age - 50) +
        I(dx_logpsa - 2) + 
        tx_grade_group + 
        tx_tstage + 
        sms + 
        lni + 
        tt_bcr + 
        st_ht + 
        st_rt + 
        slope_logpsa_st:salvage_indicator + 
        st_ht:st_rt
    formula_alpha_mets <- ~ -1 +
        slope_logpsa:as.numeric(salvage_indicator==0)
} else if (ff==3) {
    formula_gamma_mets <- ~ 1 +
        I(bcr_age - 50) +
        I(dx_logpsa - 2) + 
        tx_grade_group + 
        tx_tstage + 
        sms + 
        lni + 
        tt_bcr + 
        st_ht + 
        st_rt + 
        value_logpsa_st:salvage_indicator +
        slope_logpsa_st:salvage_indicator + 
        st_ht:st_rt
    formula_alpha_mets <- ~ -1 +
        value_logpsa:as.numeric(salvage_indicator==0) + 
        slope_logpsa:as.numeric(salvage_indicator==0)
}

formula_gamma_death <- ~ 1 +
    I(bcr_age - 50) + 
    cci +
    st_ht + 
    st_rt + 
    st_ht:st_rt

# Spline design
bs_knots <- mod_jm$control$knots[[1]]
bs_ord <- mod_jm$control$Bsplines_degree + 1

# Model expects all values to be finite, use this value to indicate treatment 
# never starts (be sure to update if time scale changes)
infty <- 100

##############################
# Extract sample of MCMC parameters and random effects, conditional on observed data
##############################
get_mcmc_sample <- function(mod_jm, data_long, data_surv, time_now, n_samples=500) {
    pred_long <- predict(mod_jm, 
                         newdata = list(newdataL = data_long, newdataE = data_surv),
                         process = 'longitudinal',
                         return_mcmc = TRUE,
                         return_newdata = TRUE,
                         times = time_now+1,
                         n_samples = n_samples,
                         n_mcmc = n_samples)
    mcmc <- attr(pred_long, 'mcmc')
    mcmc$b <- t(mcmc$b[1,,])
    mcmc$betas <- mcmc$betas[[1]]
    mcmc$D <- NULL
    
    # Split gamma coefficients into mets and death related coefficients
    gammas_mets <- mcmc$gammas[,str_detect(colnames(mcmc$gammas), 'mets')]
    colnames(gammas_mets) <- str_remove(colnames(gammas_mets), ':*dummy\\(.*?\\):*')
    mcmc$gammas_mets <- gammas_mets
    gammas_death <- mcmc$gammas[,str_detect(colnames(mcmc$gammas), 'death')]
    colnames(gammas_death) <- str_remove(colnames(gammas_death), ':*dummy\\(.*?\\):*')
    mcmc$gammas_death <- gammas_death
    mcmc$gammas <- NULL
    
    # Extract constant, linear and quadratic fixed effects coefficients
    i1 <- str_detect(colnames(mcmc$betas), 'psa_time(?!\\^2)')
    i2 <- str_detect(colnames(mcmc$betas), 'psa_time\\^2')
    mcmc$betas_0 <- mcmc$betas[, !i1 & !i2, drop=FALSE]
    mcmc$betas_1 <- mcmc$betas[, i1, drop=FALSE]
    mcmc$betas_2 <- mcmc$betas[, i2, drop=FALSE]
    
    mcmc
}
# Get one iteration of the sample
slice_mcmc <- function(mcmc, i) {
    lapply(mcmc, \(x) if(is.null(dim(x))) x else x[i,,drop=FALSE])
}

##############################
# Compute eta (or a transformation of it) at a given set of times for one
# value of b and the sampled model parameters
##############################
psa_terms <- function(data_id, mcmc) {
    a <- c(model.matrix(formula_beta_0, data_id) %*% t(mcmc$betas_0) + mcmc$b[,1])
    b <- c(model.matrix(formula_beta_1, data_id) %*% t(mcmc$betas_1) + mcmc$b[,2])
    c <- c(model.matrix(formula_beta_2, data_id) %*% t(mcmc$betas_2) + mcmc$b[,3])
    return(c(a=a, b=b, c=c))
}
psa_values <- function(terms, times, ff='value') {
    if (ff=='value') {
        return(terms['a'] + terms['b']*times + terms['c']*times^2)
    } else if (ff=='slope') {
        return(terms['b'] + 2*terms['c']*times)
    }
}

##############################
# Baseline hazard functions
##############################
bh <- function(t, mcmc, risk = 'mets', exponentiate = TRUE) {
    idx <- if(risk == 'mets') 12:22 else 1:11
    bs_gammas <- mcmc$bs_gammas[idx]
    W0 <- splineDesign(knots = bs_knots,
                       x = t,
                       ord = bs_ord,
                       outer.ok = TRUE)
    log_bh <- c(W0 %*% bs_gammas) - c(mcmc$W_std_gammas) - c(mcmc$Wlong_std_alphas)
    if(exponentiate) exp(log_bh) else log_bh
}

##############################
# Hazard functions
##############################
hazard <- function(times, 
                   data_id, 
                   mcmc, 
                   tt_salvage_tx_new, 
                   terms,
                   type=c('both', 'mets', 'death'), 
                   use_ht=1, 
                   use_rt=1) {
    type <- match.arg(type)
    
    
    dset <- cbind(data_id, data.frame(psa_time=times))
    dset <- dset |> mutate(tt_salvage_tx = tt_salvage_tx_new)
    dset <- dset |> mutate(salvage_indicator = 1*(psa_time >= tt_salvage_tx))
    dset <- dset |> mutate(st_ht = 1*salvage_indicator*use_ht, st_rt = 1*salvage_indicator*use_rt)
    dset$value_logpsa_st <- as.vector(psa_values(terms, tt_salvage_tx_new, ff='value'))
    dset$slope_logpsa_st <- as.vector(psa_values(terms, tt_salvage_tx_new, ff='slope'))
    dset$value_logpsa <- as.vector(psa_values(terms, times, ff='value'))
    dset$slope_logpsa <- as.vector(psa_values(terms, times, ff='slope'))
    
    if(type=='mets') {
        log_hr <- mcmc$gammas_mets %*% t(model.matrix(formula_gamma_mets, dset)[,-1]) +
            mcmc$alphas[,-1] %*% t(model.matrix(formula_alpha_mets, dset))
        return(exp(bh(times, mcmc, risk='mets', exponentiate=FALSE) + log_hr))
    } else if (type=='death') {
        log_hr <- mcmc$gammas_death %*% t(model.matrix(formula_gamma_death, dset)[,-1]) 
        return(exp(bh(times, mcmc, risk='death', exponentiate=FALSE) + log_hr))
    } else if(type=='both') {
        log_hr_mets <- mcmc$gammas_mets %*% t(model.matrix(formula_gamma_mets, dset)[,-1]) +
            mcmc$alphas[,-1] %*% t(model.matrix(formula_alpha_mets, dset))
        log_hr_death <- mcmc$gammas_death %*% t(model.matrix(formula_gamma_death, dset)[,-1]) 
        return(exp(bh(times, mcmc, risk='mets', exponentiate=FALSE) + log_hr_mets) + 
                   exp(bh(times, mcmc, risk='death', exponentiate=FALSE) + log_hr_death))
    }
}

# Integrate a function with a discontinuity. This is faster when integrating the
# hazard functions which have a jump discontinuity at time of treatment.
integrate_discont <- function(f, lower, upper, discont=NULL) {
    if (!is.null(discont) & lower < discont & discont < upper) {
        return(integrate(f, lower, discont)$value + integrate(f, discont, upper)$value)
    } else {
        return(integrate(f, lower, upper)$value)
    }
}

# For a numeric vector (t_1,...t_n), evaluate \int_{t_1}^{t_k} f(t)dt for each k, 
# then apply cumsum to get \int_0^{t_k} f(t)dt for each k. Slightly faster than 
# evaluating each integral separately.
integrate_along <- function(f, ts, discont=NULL) {
    int <- function(i) {
        if (i==1) return(0)
        t0 <- ts[i-1]
        t1 <- ts[i]
        integrate_discont(f, t0, t1, discont)
    }
    result <- map_dbl(seq_along(ts), int)
    cumsum(result)
}

# Compute risks of mets and death for a given treatment time and one set of
# model parameters and random effects
risk_calculation <- function(data_id,
                             mcmc,
                             time_now,
                             time_end=time_now+5,
                             tt_salvage_tx_new=infty,
                             length_out = 10,
                             length_grid = 10) {
    # Risks of mets and death will be computed at a sequence of points
    t_out <- seq(time_now, time_end, length.out=length_out)
    
    # Calculate terms of PSA trajectory
    terms <- psa_terms(data_id, mcmc)
    
    # Calculate overall survival (no mets and no death) at a grid of points. In 
    # subsequent integrals, we linearly interpolate between these values to improve
    # speed. Set the length of points `length_grid` so that the approximation error is small.
    t_grid <- seq(time_now, time_end, length.out=length_grid)
    haz_integrand <- function(s) hazard(s, data_id, mcmc, tt_salvage_tx_new, terms, type='both')
    cum_haz <- integrate_along(haz_integrand, t_grid, discont=tt_salvage_tx_new)
    surv_fn <- function(t) approx(t_grid, exp(-cum_haz), xout = t, method = 'linear', yleft = 1)$y
    
    # Probability of metastasis
    pr_mets_integrand <- function(s){
        surv_fn(s) * hazard(s, data_id, mcmc, tt_salvage_tx_new, terms, type='mets')
    }
    pr_mets <- integrate_along(pr_mets_integrand, t_out, discont=tt_salvage_tx_new)
    
    # Probability of death
    pr_death_integrand <- function(s){
        surv_fn(s) * hazard(s, data_id, mcmc, tt_salvage_tx_new, terms, type='death')
    }
    pr_death <- integrate_along(pr_death_integrand, t_out, discont=tt_salvage_tx_new)
    
    # Calculate restricted time spent in healthy+alive (ha) state and healthy+alive+untreated (hau) state
    rst_ha_integrand <- function(t) {
        (1 - approx(t_out, pr_mets, method = 'linear', xout = t)$y) * 
            (1 - approx(t_out, pr_death, method = 'linear', xout = t)$y)
    }
    rst_ha <- integrate(rst_ha_integrand, time_now, time_end)$value
    rst_hau_integrand <- function(t) {
        (1 - approx(t_out, pr_mets, method = 'linear', xout = t)$y) * 
            (1 - approx(t_out, pr_death, method = 'linear', xout = t)$y) *
            (t < tt_salvage_tx_new)
    }
    rst_hau <- integrate(rst_hau_integrand, time_now, time_end)$value
    
    return(list(risks=data.frame(time=t_out, pr_mets=pr_mets, pr_death=pr_death),
                rst_ha=rst_ha,
                rst_hau=rst_hau))
}

##############################
# Run simulations for one patient and a set of strategies
##############################
run_simulations <- function(data_id, 
                            data_psa, 
                            strategies,
                            time_now=max(data_psa$psa_time),
                            time_horizon=5,
                            n_iter=500,
                            parallel=TRUE) {
    # Create survival data in format expected by JMbayes2
    data_surv <- data_id |> mutate(
        start = 0,
        stop = time_now,
        event = 0,
        salvage_indicator = 0,
        st_ht = 0,
        st_rt = 0,
    )
    data_surv <- crossing(data_surv, tibble(cr=c('mets', 'death')))
    
    # Create longitudinal data in format expected by JMbayes2
    data_long <- crossing(data_id,  data_psa)
    data_long <- data_long |> mutate(logpsa = log(psa + 1))
    
    # Sample distribution of model parameters and random effects given data
    mcmc <- get_mcmc_sample(mod_jm, 
                            data_long, 
                            data_surv, 
                            time_now=time_now,
                            n_samples=n_iter)
    
    # Create simulated data. The indices are:
    #   i: MC simulation iteration
    #   j: patient history
    #   k: treatment strategy
    sim_data <- crossing(i = 1:n_iter,
                         k = strategies$k)
    
    # Attach MCMC sample (and random effects) for each iteration
    sim_data <- sim_data |> mutate(mcmc = map(i, \(i) slice_mcmc(mcmc, i)))
    
    # Calculate salvage therapy times for each strategy. For dynamic strategy, PSA will be measured
    # every `1/psa_frequency` years. Simulated PSA values are eta(t) + rnorm(0, sigma). Treatment 
    # time is the time of first measurement greater than `psa_threshold`
    f <- function(k, mcmc, psa_frequency=2) {
        strategy <- strategies$strategy[[k]]
        psa_threshold <- strategies$psa_threshold[[k]]
        
        if (strategy=='never') return(infty)
        else if (strategy=='now') return(time_now)
        else if (strategy=='dynamic') {
            test_times <- seq(time_now+1/psa_frequency, time_now+time_horizon, by=1/psa_frequency)
            log_psa_threshold <- log(1+psa_threshold)
            terms <- psa_terms(data_id, mcmc)
            eta <- psa_values(terms, test_times)
            y <- eta + rnorm(length(eta), 0, mcmc$sigmas)
            tt_salvage_tx <- if(all(y < log_psa_threshold)) infty else test_times[min(which(y >= log_psa_threshold))]
            return(tt_salvage_tx)
        }
    }
    sim_data$tt_salvage_tx <- map2_dbl(sim_data$k, sim_data$mcmc, f)
    
    # Compute risk of metastasis and death for each iteration
    f <- function(i, k, mcmc, tt_salvage_tx) risk_calculation(data_id = data_id,
                                                                 mcmc = mcmc,
                                                                 time_now = time_now,
                                                                 time_end = time_now + time_horizon,
                                                                 tt_salvage_tx_new = tt_salvage_tx)
    map_fn <- if (parallel) future_pmap else pmap
    sim_data$out <- map_fn(sim_data, f)
    sim_data <- sim_data |> unnest_wider(out)
    
    # Compute mean PSA values (eta) over follow-up time for each iteration
    f <- function(k, mcmc, d_out=0.1, ...) {
        time <- seq(0, time_now+time_horizon, by=d_out)
        logpsa <- psa_values(psa_terms(data_id, mcmc), time, ff='value')
        return(data.frame(time=time, logpsa=logpsa))
    }
    sim_data$psa <- pmap(sim_data, f)
    
    # Aggregate across iterations
    
    # Calculate mean risk across iterations
    rset <- sim_data |> unnest(risks)
    rset <- rset |> group_by(k, time)
    rset <- rset |> summarize(pr_mets = mean(pr_mets),
                              pr_death = mean(pr_death))
    
    # Calculate quantiles of PSA trajectories
    pset <- sim_data |> filter(k==1) |> unnest(psa)
    pset <- pset |> group_by(time)
    pset <- pset |> summarize(q25=quantile(logpsa, 0.25),
                              q50=quantile(logpsa, 0.50),
                              q75=quantile(logpsa, 0.75))
    
    # Calculate probability of starting treatment
    f <- function(t) {
        u <- c(0, sort(unique(t)))
        p <- map_dbl(u, \(x) sum(t <= x)) / length(t)
        return(data.frame(time = u, pr_tx = p))
    }
    sset <- sim_data |> group_by(k)
    sset <- sset |> reframe(tx = f(tt_salvage_tx)) |> unnest(tx)
    
    # Calculate restricted mean time in states
    tset <- sim_data |> group_by(k)
    tset <- tset |> summarize(rst_ha = mean(rst_ha),
                              rst_hau = mean(rst_hau))
    
    list(rset=rset,
         pset=pset,
         sset=sset,
         tset=tset)
}

data_id <- tibble(
    ptid = '1',
    bcr_age = 75,
    cci = cci_fct('2+'),
    dx_logpsa = log(1+4),
    tx_grade_group = tx_grade_group_fct('4+3'),
    tx_tstage = tx_tstage_fct('T3a'),
    sms = sms_lni_fct('Negative'),
    lni = sms_lni_fct('Negative'),
    tt_bcr = 2,
    tt_salvage_tx = infty,
    value_logpsa_st = 0,
    slope_logpsa_st = 0
)

data_psa_1 <- tribble(~psa_time, ~psa,
                    0.00, 0.20,
                    0.25, 0.20)
data_psa_2 <- tribble(~psa_time, ~psa,
                      0.00, 0.20,
                      0.25, 0.20,
                      0.50, 0.25,
                      0.75, 0.25,
                      1.00, 0.25,
                      1.25, 0.25)


# Treatment strategies to consider
strategies <- tribble(
    ~k, ~strategy, ~psa_threshold,
    1,  'never',   NA,
    2,  'now',     NA,
    3,  'dynamic', 0.5,
    4,  'dynamic', 1.0,
    5,  'dynamic', 2.0
)
strategies$k <- as.integer(strategies$k)

results_1 <- run_simulations(data_id, data_psa_1, strategies)
results_2 <- run_simulations(data_id, data_psa_2, strategies)

##############################
# Plot results
##############################
plot_results <- function(results, time_now) {
    theme_set(theme_classic(base_size = 8))
    
    time_end <- time_now + 5
    
    # PSA data
    pset <- results$pset
    
    # Treatment data
    sset <- results$sset |> filter(k>2)
    sset <- sset |> left_join(strategies)
    sset <- sset |> mutate(label = str_glue('{scales::number(psa_threshold, 0.1)} ng/mL'))
    sset <- sset |> filter(time < time_end)
    sset <- rbind(sset, sset |> group_by(k) |>  slice_tail(n=1) |> mutate(time=time_end))
    
    # Risk data
    rset <- results$rset 
    rset <- rset |> left_join(strategies)
    rset <- rset |> mutate(label = case_when(strategy == 'now' ~ 'Immediate',
                                             strategy == 'never' ~ 'Never',
                                             strategy == 'dynamic' ~ str_glue('Dynamic (PSA threshold {scales::number(psa_threshold, 0.1)} ng/mL)')))
    rset <- rset |> mutate(label = fct_inorder(label))
    
    # Restricted mean times
    tset <- results$tset
    tset <- tset |> left_join(strategies)
    tset <- tset |> mutate(label = case_when(strategy=='never' ~ 'Never',
                                             strategy=='now' ~ 'Immediate',
                                             strategy=='dynamic' ~ str_glue('Threshold {scales::number(psa_threshold, 0.1)} ng/mL')))
    tset <- tset |> mutate(label = fct_reorder(label, rst_ha))
    
    # PSA trajectories
    ymax <- 3
    pset <- pset |> mutate(across(starts_with('q'), \(x) pmax(pmin(x, ymax), 0)))
    g1 <- ggplot(pset)
    g1 <- g1 + geom_ribbon(aes(x = time, ymin = q25, ymax = q75), 
                           alpha = 0.25, 
                           fill = 'lightblue', 
                           color = 'darkgray')
    g1 <- g1 + geom_line(aes(x = time, y = q50), 
                         color = 'darkblue', 
                         linewidth = 0.5)
    # g1 <- g1 + geom_point(data = patients$data_psa[[j0]], 
    #                       aes(x = psa_time, y = log(psa+1)), 
    #                       shape = 1, 
    #                       size = 1.5)
    g1 <- g1 + geom_vline(aes(xintercept = time_now), linetype = 'dashed')
    g1 <- g1 + geom_hline(aes(yintercept = log(1+0.5)), linetype = 'dotted')
    g1 <- g1 + geom_hline(aes(yintercept = log(1+1.0)), linetype = 'dotted')
    g1 <- g1 + geom_hline(aes(yintercept = log(1+2.0)), linetype = 'dotted')
    g1 <- g1 + annotate('text', x = time_now+0.1, y = 2.5, 
                        label = str_glue('Current time:\n{time_now} years'),
                        hjust = 'left', size = 2)
    g1 <- g1 + annotate('text', x = time_end, y = log(1.5)+0.02, vjust = 'bottom', hjust = 'right', label = 'PSA 0.5ng/mL', size = 2)
    g1 <- g1 + annotate('text', x = time_end, y = log(2.0)+0.02, vjust = 'bottom', hjust = 'right', label = 'PSA 1.0ng/mL', size = 2)
    g1 <- g1 + annotate('text', x = time_end, y = log(3.0)+0.02, vjust = 'bottom', hjust = 'right', label = 'PSA 2.0ng/mL', size = 2)
    g1 <- g1 + scale_x_continuous(expand = c(0,0), limits = c(0, time_end), breaks = seq(0, time_end, by=1))
    g1 <- g1 + scale_y_continuous(expand = c(0,0), limits = c(0,ymax))
    g1 <- g1 + labs(x = 'Years after BCR', y = 'log(PSA + 1)')
    g1
    
    # Probability of starting treatment
    
    g2 <- ggplot(sset)
    g2 <- g2 + geom_step(aes(x = time, y = pr_tx, linetype=label),
                         color = '#009E73', 
                         linewidth = 0.5)
    g2 <- g2 + geom_vline(aes(xintercept = time_now), 
                          linetype = 'dashed')
    g2 <- g2 + scale_x_continuous(name = 'Years after BCR',
                                  expand = c(0,0), 
                                  limits = c(0, time_end), 
                                  breaks = seq(0, time_end, by=1))
    g2 <- g2 + scale_y_continuous(name = 'Probability of Treatment',
                                  expand = c(0,0), 
                                  limits = c(0, 1.05), 
                                  labels = scales::percent)
    g2 <- g2 + scale_linetype_manual(name = 'PSA threshold',
                                     values = c('solid', 'dashed', 'dotted'))
    g2 <- g2 + theme(legend.position = c(0.98, 0.02),
                     legend.justification = c('right', 'bottom'),
                     legend.text = element_text(margin = margin(t=0, b=0)),
                     legend.key.height = unit(0.1, 'cm'),
                     legend.spacing.y = unit(0.1, 'cm'),
                     legend.background = element_rect(color = 'black'))
    g2
    
    g3 <- ggplot(rset)
    g3 <- g3 + geom_line(aes(x=time, y=pr_mets, color=label, linetype=label),
                         linewidth=0.5)
    g3 <- g3 + geom_vline(aes(xintercept = time_now), linetype = 'dashed')
    g3 <- g3 + scale_x_continuous(expand = c(0,0), limits = c(0, time_end), breaks = seq(0, time_end, by=1))
    g3 <- g3 + scale_y_continuous(expand = c(0,0), limits = c(0,0.51), breaks=seq(0,0.5, by=0.1), labels = scales::percent)
    g3 <- g3 + scale_color_manual(values = c('#000000', '#009E73', '#E69F00', '#E69F00',  '#E69F00'))
    g3 <- g3 + scale_linetype_manual(values = c('solid', 'solid', 'solid', 'dashed', 'dotted'))
    g3 <- g3 + labs(x = 'Years after BCR', y = 'Probability of Metastasis',
                    color = 'Treatment strategy', linetype = 'Treatment strategy')
    g3 <- g3 + theme(legend.position = c(0.02, 0.98),
                     legend.justification = c('left', 'top'),
                     legend.text = element_text(margin = margin(t=0, b=0)),
                     legend.key.height = unit(0.1, 'cm'),
                     legend.spacing.y = unit(0.1, 'cm'),
                     legend.background = element_rect(color = 'black'))
    g3
    
    g4 <- ggplot(rset)
    g4 <- g4 + geom_line(aes(x=time, y=pr_death, color=label, linetype=label),
                         linewidth=0.5)
    g4 <- g4 + geom_vline(aes(xintercept = time_now), linetype = 'dashed')
    g4 <- g4 + scale_x_continuous(expand = c(0,0), limits = c(0, time_end), breaks = seq(0, time_end, by=1))
    g4 <- g4 + scale_y_continuous(expand = c(0,0), limits = c(0,0.51), breaks=seq(0,0.5, by=0.1), labels = scales::percent)
    g4 <- g4 + scale_color_manual(values = c('#000000', '#009E73', '#E69F00', '#E69F00',  '#E69F00'))
    g4 <- g4 + scale_linetype_manual(values = c('solid', 'solid', 'solid', 'dashed', 'dotted'))
    g4 <- g4 + labs(x = 'Years after BCR', y = 'Probability of Death',
                    color = 'Treatment strategy', linetype = 'Treatment strategy')
    g4 <- g4 + theme(legend.position = 'none')
    g4
    
    g <- cowplot::plot_grid(g1, g2, g3, g4, ncol=2, labels=c('A', 'B', 'C', 'D'), label_size = 10, align = 'hv')
    g
    
    # Restricted mean times
    tset <- tset |> mutate(rst_ha = rst_ha-rst_hau)
    tset <- tset |> pivot_longer(cols = starts_with('rst'), names_prefix = 'rst_')
    g5 <- ggplot(tset)
    g5 <- g5 + geom_col(aes(x = 5, y = label), 
                        fill = NA, 
                        color = 'black', 
                        position = 'dodge', 
                        linewidth = 0.25, 
                        width = 0.75)
    g5 <- g5 + geom_col(aes(x = value, y = label, fill = name), 
                        position = 'stack', 
                        color = 'black', 
                        linewidth = 0.25, 
                        width = 0.75)
    g5 <- g5 + scale_fill_manual(values = c('#56B4E9', '#0072B2'),
                                 labels = c('Alive, unprogressed', 'Alive, unprogressed, untreated'))
    g5 <- g5 + labs(x = '5-year restricted mean time in state', y = 'Treatment strategy', fill = 'State')
    g5 <- g5 + scale_x_continuous(limits = c(0, 5.1), expand = c(0,0))
    g5 <- g5 + guides(fill = guide_legend(keywidth = 0.75, keyheight = 0.75))
    g5 <- g5 + theme(legend.position = 'bottom')
    g5
    
    g <- cowplot::plot_grid(g, g5, ncol=1, labels = c('', 'E'), label_size = 10, rel_heights = c(2,0.8))
    g
    
}
plot_results(results_1, time_now=0.25)
plot_results(results_2, time_now=1.25)

