# main.R script for C-LARS
# ------------------------
rm(list = ls())

source("generate_cmplx_gauss_data.R")
source("complex_lars.R")

generate_simdata <- FALSE
import_ws <- TRUE
ws_small <- TRUE

if (generate_simdata) {
  n <- 100
  p <- 150000
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


# employ CLARS
# ------------------------------
tictoc::tic()
res <- complex_lars(y = data$y,
                    x = data$X,
                    has_intercept = FALSE,
                    standardize = TRUE,
                    desired_model_size = 5,
                    tol = .Machine$double.eps)
tictoc::toc()


# Test
# -----------------------------
cat("Proposed index set by CLARS: ", res$index_actives, "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
    intersect(res$index_actives, drop(data$support)), "\n")

# Test other output parameters
# -----------------------------
res$beta_zm
res$y_hat
res$ssr
res$r2
