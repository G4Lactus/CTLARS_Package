rm(list = ls())

source("generate_cmplx_gauss_data.R")
source("obj_complex_lars.R")

generate_simdata <- FALSE
import_ws <- TRUE
ws_small <- TRUE

if (generate_simdata) {
  n <- 100
  p <- 1500
  beta_cardinality <- 5
  set.seed(42)

  # generate some data
  data <- generate_cmplx_gauss_data(
    n_rows = n,
    p_cols = p,
    beta_cardinality = beta_cardinality,
    scenario = "standard",
    beta_type = "complex", # "complex", "real", "imaginary",
    beta_val = "random" # "constant", "constant"
  )
} else {
  library(R.matlab)
  if (import_ws) {
    if (ws_small) {
      data <- R.matlab::readMat("WS_CLARS_small.mat")
    } else {
      data <- R.matlab::readMat("WS_CLARS_medium.mat")
    }
  }
}


tictoc::tic()
clars_obj <- clars$new(data$X, data$y, has_intercept = FALSE)
clars_obj$execute_clars_step(desired_model_size = 5)
tictoc::toc()


# Test
# -----------------------------
cat("Proposed index set by CLARS: ", clars_obj$get_active_set(), "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
    intersect(clars_obj$get_active_set(), drop(data$support)), "\n")
