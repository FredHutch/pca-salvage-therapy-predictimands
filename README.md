# pca-salvage-therapy-predictimands

Code supporting "Predictimands for Salvage Therapy in Biochemically Recurrent Prostate Cancer Using Joint Longitudinal and Failure Time Models"

## Scripts

**`01-fit-models.R`** — Fits the joint longitudinal-survival model to the MSK
clinical dataset. Runs linear mixed-effects and competing-risks Cox submodels,
then combines them into a Bayesian joint model using JMbayes2.

**`02-dynamic-predictions.R`** — Computes dynamic predictions and predictimands
from the fitted joint model estimated in the previous step. For a patient
observed up to a landmark time, estimates the probability of metastasis under
three salvage therapy strategies (immediate, never, dynamic threshold-based) by
integrating over posterior draws of random effects and model parameters.

**`03-simulation-study.R`** — Validates the predictimand estimation approach
via simulation. Step 1 generates synthetic longitudinal and competing-risks
data and fits a joint model to each replicate. Step 2 computes true and
estimated predictimands for a cohort of simulated patients, using the procedure
outlined in the paper.


