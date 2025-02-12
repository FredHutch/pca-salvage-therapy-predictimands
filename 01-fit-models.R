###################################
# 01-fit-models.R
#
# Fit a joint model for post-BCR PSA evolution and competing risks of 
# metastatic progression and other-cause death
#
# PSA: Prostate-specific antigen
# ST: Salvage therapy
# RP: Radical prostatectomy
# BCR: Biochemical recurrence
###################################

library(JMbayes2)
library(tidyverse)

#' Starting data
#' The starting data consists of `data_long` with longitudinal PSA measurements
#' and `data_surv` with competing risks data. The formats of these files is below
#'
#' `data_long`
#'  Rows: 1 row per patient and PSA observation
#'  Columns:
#'      ptid: Patient identifier
#'      psa_time: Time of PSA lab, measured in years since BCR
#'      logpsa: log(PSA+1)
#'      tx_grade_group: Pathology grade (<=6, 3+4, 4+3, 8, 9-10)
#'      tx_tstage: Pathology T-stage (T0-T2, T3a, T3b, T4)
#'      lni: Lymph node involvement
#'      sms: Surgical margins
#'      dx_logpsa: log(PSA+1) at time of RP
#'      tt_bcr: Time from RP to BCR
#'       
#' `data_surv`
#' Rows: Up to four rows per patient, representing each competing risk and before/after ST
#' Columns:
#'      ptid, tx_grade_group, tx_tstage, lni, sms, dx_logpsa, tt_bcr: Same as above
#'      bcr_age: Age at BCR
#'      cci: Charlson comorbidity index
#'      cr: Row identifier for competing risks ('mets' or 'death')
#'      salvage_indicator: Time-varying indicator of salvage therapy
#'      st_ht: Time-varying indicator of hormone therapy
#'      st_rt: Time-varying indicator of radiation therapy
#'
#'       

###################################
# Fit longitudinal submodel
###################################
mod_long <- lme(data = data_long,
                control = nlme::lmeControl(opt = 'optim'),
                method = 'ML',
                fixed = logpsa ~ 1 +
                    # Intercept
                    tx_grade_group +
                    lni +
                    I(dx_logpsa - 2) +
                    tt_bcr +
                    # Linear time
                    psa_time + 
                    psa_time:tx_grade_group +
                    psa_time:lni +
                    psa_time:tt_bcr + 
                    # Quadratic time
                    I(psa_time^2),
                random = list(ptid = pdDiag(form = ~ 1 + psa_time + I(psa_time^2))))
summary(mod_long)


###################################
# Attach fitted values and slopes of PSA
###################################
# Extract random effects for each patient
dset <- data_surv |> group_by(ptid) |> slice_head(n = 1)
random_effects <- ranef(mod_long)
names(random_effects) <- c('b0', 'b1', 'b2')
random_effects <- as_tibble(random_effects, rownames = 'ptid')
dset <- dset |> left_join(random_effects)

# Calculate intercept, slope and quadratic time terms (incl both fixed and random effects in each)
beta <- mod_long$coefficients$fixed
m0 <- model.matrix(~ tx_grade_group + lni + I(dx_logpsa - 2) + tt_bcr, dset)
dset$c0 <- as.vector(m0 %*% beta[colnames(m0)]) + dset$b0
m1 <- model.matrix(~ 1 + tx_grade_group + lni + tt_bcr, dset)
beta1 <- beta[ifelse(colnames(m1)=='(Intercept)', 'psa_time', paste0(colnames(m1), ':psa_time'))]
dset$c1 <- as.vector(m1 %*% beta1) + dset$b1
dset$c2 <- as.vector(beta['I(psa_time^2)']) + dset$b2

# Calculate PSA value and slope
dset <- dset |> mutate(
    value_logpsa_st = ifelse(tt_salvage_tx > 999, 0, c0 + c1*tt_salvage_tx + c2*tt_salvage_tx^2),
    slope_logpsa_st = ifelse(tt_salvage_tx > 999, 0, c1 + 2*c2*tt_salvage_tx),
    area_logpsa_st = ifelse(tt_salvage_tx > 999, 0, c0 + c1*tt_salvage_tx/2 + c2*tt_salvage_tx^2/3),
)
data_surv <- left_join(data_surv, select(dset, ptid, value_logpsa_st, slope_logpsa_st, area_logpsa_st))


###################################
# Fit hazard submodel, and then joint model, for each selected functional form
###################################
n_iter <- 50000
# ff <- '1'
# ff <- '2'
ff <- '3'

dummy <- function (f, lvl) as.numeric(f == lvl)
if (ff == '1') {
    functional_forms <- ~ dummy(cr, 'death'):value(logpsa) +
        dummy(cr, 'mets'):dummy(salvage_indicator, 0):value(logpsa)
    priors <- list(Tau_alphas = list(diag(50000, 1), diag(1)),
                   gamma_prior_D_sds = FALSE)
    surv_formula <- Surv(start, stop, event) ~ strata(cr) +
        # Metastasis risk
        dummy(cr, 'mets'):I(bcr_age - 50) + 
        dummy(cr, 'mets'):I(dx_logpsa - 2) +
        dummy(cr, 'mets'):tx_grade_group +
        dummy(cr, 'mets'):tx_tstage +
        dummy(cr, 'mets'):sms +
        dummy(cr, 'mets'):lni +
        dummy(cr, 'mets'):tt_bcr +
        dummy(cr, 'mets'):value_logpsa_st:salvage_indicator +
        dummy(cr, 'mets'):st_ht + 
        dummy(cr, 'mets'):st_rt + 
        dummy(cr, 'mets'):st_ht:st_rt + 
        # Other-cause death risk
        dummy(cr, 'death'):I(bcr_age - 50) + 
        dummy(cr, 'death'):cci + 
        dummy(cr, 'death'):st_ht + 
        dummy(cr, 'death'):st_rt + 
        dummy(cr, 'death'):st_ht:st_rt
} 
if (ff == '2') {
    functional_forms <- ~ dummy(cr, 'death'):value(logpsa) +
        dummy(cr, 'mets'):dummy(salvage_indicator, 0):slope(logpsa)
    priors <- list(Tau_alphas = list(diag(50000, 1), diag(1)),
                   gamma_prior_D_sds = FALSE)
    surv_formula <- Surv(start, stop, event) ~ strata(cr) +
        # Metastasis risk
        dummy(cr, 'mets'):I(bcr_age - 50) + 
        dummy(cr, 'mets'):I(dx_logpsa - 2) +
        dummy(cr, 'mets'):tx_grade_group +
        dummy(cr, 'mets'):tx_tstage +
        dummy(cr, 'mets'):sms +
        dummy(cr, 'mets'):lni +
        dummy(cr, 'mets'):tt_bcr +
        dummy(cr, 'mets'):slope_logpsa_st:salvage_indicator +
        dummy(cr, 'mets'):st_ht + 
        dummy(cr, 'mets'):st_rt + 
        dummy(cr, 'mets'):st_ht:st_rt + 
        # Other-cause death risk
        dummy(cr, 'death'):I(bcr_age - 50) + 
        dummy(cr, 'death'):cci + 
        dummy(cr, 'death'):st_ht + 
        dummy(cr, 'death'):st_rt + 
        dummy(cr, 'death'):st_ht:st_rt
} 
if (ff == '3') {
    functional_forms <- ~ dummy(cr, 'death'):value(logpsa) +
        dummy(cr, 'mets'):dummy(salvage_indicator, 0):value(logpsa) +
        dummy(cr, 'mets'):dummy(salvage_indicator, 0):slope(logpsa)
    priors <- list(Tau_alphas = list(diag(50000, 1), diag(1), diag(1)),
                   gamma_prior_D_sds = FALSE)
    surv_formula <- Surv(start, stop, event) ~ strata(cr) +
        # Metastasis risk
        dummy(cr, 'mets'):I(bcr_age - 50) + 
        dummy(cr, 'mets'):I(dx_logpsa - 2) +
        dummy(cr, 'mets'):tx_grade_group +
        dummy(cr, 'mets'):tx_tstage +
        dummy(cr, 'mets'):sms +
        dummy(cr, 'mets'):lni +
        dummy(cr, 'mets'):tt_bcr +
        dummy(cr, 'mets'):value_logpsa_st:salvage_indicator +
        dummy(cr, 'mets'):slope_logpsa_st:salvage_indicator +
        dummy(cr, 'mets'):st_ht + 
        dummy(cr, 'mets'):st_rt + 
        dummy(cr, 'mets'):st_ht:st_rt + 
        # Other-cause death risk
        dummy(cr, 'death'):I(bcr_age - 50) + 
        dummy(cr, 'death'):cci + 
        dummy(cr, 'death'):st_ht + 
        dummy(cr, 'death'):st_rt + 
        dummy(cr, 'death'):st_ht:st_rt
}

mod_surv <- coxph(data = data_surv,
                  formula = surv_formula)

mod_jm <- jm(mod_surv, 
             mod_long, 
             data_Surv = data_surv, 
             data_Long = data_long,
             time_var = 'psa_time',  
             parallel = 'multicore',
             n_iter = n_iter,
             n_burnin = as.integer(n_iter/10),
             functional_forms = functional_forms,
             priors = priors)
summary(mod_jm)

# Save output
save(mod_surv, file = str_glue(str_glue('output/model-surv-ff{ff}.RData')))
save(mod_jm, file = str_glue('output/model-jm-ff{ff}.RData'))
