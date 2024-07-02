rm(list = ls())

source("generate_cmplx_gauss_data.R")
source("obj_complex_tlars.R")

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


# append dummies
p <- ncol(data$X)
n <- nrow(data$X)
num_dummies <- 1 * p
data$X <- cbind(data$X,
  matrix(complex(real = stats::rnorm(n * p),
                 imag = stats::rnorm(n * num_dummies)),
         nrow = n, ncol = num_dummies)
)


cat("Start Complex Terminating Lars\n")
tictoc::tic()
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = num_dummies)
tictoc::toc()
tictoc::tic()
ctlars_obj$execute_clars_step(t_stop = 3)
tictoc::toc()


# Test
# -----------------------------
cat("Proposed index set by CLARS: ", ctlars_obj$get_active_set(), "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
    intersect(ctlars_obj$get_active_set(), drop(data$support)), "\n")
cat("Dummies in proposal: ", ctlars_obj$get_active_dummies(), "\n")


# Iterative results
# -----------------------------
tictoc::tic()
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = num_dummies)
tictoc::toc()
tictoc::tic()
for (t in 1:10) {
  ctlars_obj$execute_clars_step(t_stop = t)
}
tictoc::toc()
cat("Proposed index set by CLARS: ", ctlars_obj$get_active_set(), "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
    intersect(ctlars_obj$get_active_set(), drop(data$support)), "\n")
cat("Dummies in proposal: ", ctlars_obj$get_active_dummies(), "\n")



# Simulate T-Rex random experiments
# -----------------------------------
set.seed(42)
num_rnd_exps <- 20
t_stop <- 3
num_dummies <- 1 * p
data <- generate_cmplx_gauss_data(
  n_rows = n,
  p_cols = p,
  beta_cardinality = beta_cardinality,
  scenario = "standard",
  beta_type = "complex", # "complex", "real", "imaginary",
  beta_val = "random" # "constant", "constant"
)
experiment_data <- lapply(seq(num_rnd_exps), function(experiment) {
  data_exp <- data
  data_exp$X <- cbind(data_exp$X,
    matrix(
      complex(real = stats::rnorm(n * p),
              imag = stats::rnorm(n * num_dummies)),
      nrow = n, ncol = num_dummies
    )
  )

  return(data_exp)
})

# Select dummies until stopping crition: t_stop is matched
for (t_stop_x in seq(t_stop)) {
  if (t_stop_x == 1) {
    ctlars_objs <- lapply(seq(num_rnd_exps), function(experiment) {
      ctlars_obj_x <- ctlars$new(experiment_data[[experiment]]$X,
                                 experiment_data[[experiment]]$y,
                                 has_intercept = FALSE,
                                 num_dummies = num_dummies)
      return(ctlars_obj_x)
    })
  }


  lapply(ctlars_objs, function(ctlars_x) {
    ctlars_x$execute_clars_step(t_stop = t_stop_x)
  })
  cat("-------------------\n")
}

# check active set
ctlars_objs[[1]]$get_active_set()
ctlars_objs[[2]]$get_active_set()
ctlars_objs[[3]]$get_active_set()
ctlars_objs[[4]]$get_active_set()
ctlars_objs[[5]]$get_active_set()

remove_dummy_inds <- function(active_set, p_cols, num_dummies) {
  return(setdiff(active_set, seq(from = p_cols + 1, to = p_cols + num_dummies)))
}
# remove dummies
active_vars <- lapply(ctlars_objs, function(ctlars_x) {
  remove_dummy_inds(ctlars_x$get_active_set(), p, num_dummies)
})

mat <- matrix(0, nrow = num_rnd_exps, ncol = max(sapply(active_vars, length)))
counter <- 1
for (row_x in 1:num_rnd_exps) {
  mat[row_x, seq_along(active_vars[[row_x]])] <- active_vars[[row_x]]
}
data$support
