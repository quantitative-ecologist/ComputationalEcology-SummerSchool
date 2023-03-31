# ========================================================================

#               Multivariate models of predator behaviour                #

# ========================================================================

# Code for a multivariate mixed-models to :
# 1. partition the variance in predator behaviour
# 2. quantify behavioural correlations among predator behaviours
# -----------------------------------------------------------------------





# =======================================================================
# 1. Load libraries, and import dataset
# =======================================================================


# Detect number of cores ------------------------------------------------
options(mc.cores = parallel::detectCores())



# Load libraries --------------------------------------------------------
library(data.table)
library(brms)
library(parallel)



# Import dataset --------------------------------------------------------

# Folder path Compute Canada
folder <- file.path("/home", "maxime11", "projects", "def-monti",
                    "maxime11", "phd_project", "data")

# Import the data
data <- fread(
    file.path(folder, "FraserFrancoetal2022-data.csv"),
    select = c(
        "player_id", "avatar_id",
        "hunting_success", "game_duration",
        "speed", "space_covered_rate", "hook_start_time",
        "prey_avg_speed", "prey_avg_space_covered_rate"
    ),
    stringsAsFactors = TRUE
)

# Rename variables
setnames(data, "hook_start_time", "latency_1st_capture")
setnames(data, "speed", "pred_speed")
setnames(data, "player_id", "predator_id")

# =======================================================================
# =======================================================================





# =======================================================================
# 2. Prepare variables for the model
# =======================================================================


# Transform --------------------------------------------------------------

data[, ":=" (
  latency_1st_capture = log(latency_1st_capture + 1),
  game_duration = sqrt(game_duration)
  )
]



# Standardize the variables ----------------------------------------------

standardize <- function (x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

data[
  , c("Zgame_duration",
      "Zpred_speed",
      "Zspace_covered_rate",
      "Zlatency_1st_capture",
      "Zprey_avg_speed", 
      "Zprey_avg_space_covered_rate") :=
  lapply(.SD, standardize),
  .SDcols = c(4:9)
]

# =======================================================================
# =======================================================================





# =======================================================================
# 3. Build the multivariate models
# =======================================================================


# Model 1 ---------------------------------------------------------------

# Model for the predator's speed
speed_form1 <- bf(
  Zpred_speed ~
  1 +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()

# Model for the rate of space covered by the predator
space_form1 <- bf(
  Zspace_covered_rate ~
  1 +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()

# Model for the latency before the 1st capture
hook_form1 <- bf(
  Zlatency_1st_capture ~
  1 +
  Zgame_duration +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()



# Model 2 ---------------------------------------------------------------

# Model for the predator's speed
speed_form2 <- bf(
  Zpred_speed ~
  1 +
  Zprey_avg_speed +
  Zprey_avg_space_covered_rate +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()

# Model for the rate of space covered by the predator
space_form2 <- bf(
  Zspace_covered_rate ~
  1 +
  Zprey_avg_speed +
  Zprey_avg_space_covered_rate +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()

# Model for the latency before the 1st capture
hook_form2 <- bf(
  Zlatency_1st_capture ~
  1 +
  Zprey_avg_speed +
  Zprey_avg_space_covered_rate +
  Zgame_duration +
  (1 | a | avatar_id) +
  (1 | b | predator_id)
) + gaussian()



# priors ----------------------------------------------------------------

# Priors for model 1
priors1 <- c(
  # priors on fixed effects
  set_prior(
    "normal(0, 2)",
    class = "b",
    coef = "Zgame_duration",
    resp = "Zlatency1stcapture"
  ),
  # priors on var. parameters (brms automatically detects half-normal)
  set_prior(
    "normal(0, 1)",
    class = "sd", # applies to all variance parameters
    resp = c("Zpredspeed", "Zspacecoveredrate", "Zlatency1stcapture")
  ),
  # priors on the variance-covariance matrices
  set_prior(
    "lkj(2)",
    class = "cor",
    group = "avatar_id"
  ),
  set_prior(
    "lkj(2)",
    class = "cor",
    group = "predator_id"
  )
)

# Priors for model 2
priors2 <- c(
  priors1,
  set_prior(
    "normal(0, 2)",
    class = "b",
    coef = c("Zprey_avg_speed", "Zprey_avg_space_covered_rate"),
    resp = c("Zpredspeed", "Zspacecoveredrate", "Zlatency1stcapture",
             "Zpredspeed", "Zspacecoveredrate", "Zlatency1stcapture")
  )
)

# ======================================================================
# ======================================================================





# =======================================================================
# 4. Run the multivariate models
# =======================================================================

# Model 1 ---------------------------------------------------------------

mv_model1 <- brm(
  # The three model formulas are summed
  speed_form1 +
  space_form1 +
  hook_form1 +
  # To estimate residual correlations
  # i.e. inference about behavioural plasticity
  set_rescor(TRUE),
  # MCMC settings to obtain 1000 posterior samples
  # (iter - warmups) / thin * chains
  warmup = 500,
  iter = 2500,
  thin = 8,
  chains = 4,
  # Initialize parameter values at 0
  inits = "0",
  # Within-chain parallelization :
  # Use only if you have access to multiple computer cores
  threads = threading(12),
  # Software backend for MCMC estimation
  backend = "cmdstanr",
  seed = 123,
  prior = priors1,
  # Helps MCMC convergence
  control = list(adapt_delta = 0.95),
  # Sample priors to later compare priors vs posterior distributions
  sample_prior = TRUE,
  data = data
)

# Save the object to a specified folder path
saveRDS(mv_model1, path = "mv_model1.rds")



# Model 2 ---------------------------------------------------------------

# Model 1
mv_model2 <- brm(
  # The three model formulas are summed
  speed_form2 +
  space_form2 +
  hook_form2 +
  # To estimate residual correlations
  # i.e. inference about behavioural plasticity
  set_rescor(TRUE),
  # MCMC settings to obtain 1000 posterior samples
  # (iter - warmups) / thin * chains
  warmup = 500,
  iter = 2500,
  thin = 8,
  chains = 4,
  # Initialize parameter values at 0
  inits = "0",
  # Within-chain parallelization :
  # Use only if you have access to multiple computer cores
  threads = threading(12),
  # Software backend for MCMC estimation
  backend = "cmdstanr",
  seed = 123,
  prior = priors2,
  # Helps MCMC convergence
  control = list(adapt_delta = 0.95),
  # Sample priors to later compare priors vs posterior distributions
  sample_prior = TRUE,
  data = data
)

# Save the object to a specified folder path
saveRDS(mv_model2, path = "mv_model2.rds")

# =======================================================================
# =======================================================================