rm(list = ls()) # clear the workspace

library(ctlars)

generate_simdata <- FALSE
if (generate_simdata) {
  n <- 100
  p <- 150
  beta_cardinality <- 5
  set.seed(42)

  # generate some data
  data <- generate_ccg_data(n,
                            p,
                            mean_xr = 0,
                            sd_xr = 1,
                            mean_xi = 0,
                            sd_xi = 1,
                            beta_cardinality = beta_cardinality,
                            noise_power = 1,
                            set_snr = TRUE,
                            snr_is_linear = TRUE,
                            snr_val_linear = 2,
                            snr_val_db = 10)

} else {
  data <- ctlars::WCCGauss_data
}


# Append dummies
# ----------------------------------------
set.seed(789)
update_data <- function(data_lst, num_dummies = p) {
  data_lst$p <- ncol(data_lst$X)
  data_lst$n <- nrow(data_lst$X)

  data_lst$num_dummies <- num_dummies
  data_lst$X <- cbind(data_lst$X,
    matrix(complex(real = stats::rnorm(data_lst$n * num_dummies, sd = 1),
                   imag = stats::rnorm(data_lst$n * num_dummies, sd = 1)),
           nrow = data_lst$n, ncol = num_dummies)
  )

  return(data_lst)
}
data <- update_data(data, num_dummies = ncol(data$X))


# Run complex terminating lars algorithm
# ----------------------------------------
cat("Start Complex Terminating Lars\n")
tic_init <- tictoc::tic()
ctlars_obj <- ctlars::ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies,
                         verbose = TRUE)
toc_init <- tictoc::toc()
tic_run <- tictoc::tic()
ctlars_obj$execute_clars_step(t_stop = 3, early_stop = TRUE, use_chol = TRUE)
toc_run <- tictoc::toc()


# Test
# ----------------------------------------
cat("Proposed index set by CLARS: ", ctlars_obj$get_active_set(), "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
    intersect(ctlars_obj$get_active_set(), drop(data$support)), "\n")
cat("Dummies in proposal: ", ctlars_obj$get_active_dummies(), "\n")

## Plot
## ----------------------------------------
plot(ctlars_obj)


# Warm Start
# -------------
ctlars_obj$execute_clars_step(t_stop = 5)
plot(ctlars_obj)


## Try a second example with t_stop = 5
## ---------------------------------------
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies)

ctlars_obj$execute_clars_step(t_stop = 5, early_stop = TRUE)
plot(ctlars_obj)


# Full Lars path
# ----------------------------------------
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies)

ctlars_obj$execute_clars_step(t_stop = 1, early_stop = FALSE)
plot(ctlars_obj)

# Evaluate FDP and TPP
# ----------------------------------------
chosen <- ctlars_obj$get_active_set()
# remove dummies
chosen <- setdiff(
  chosen, ctlars_obj$get_active_dummies()
)
FDP(chosen, drop(data$support))
TPP(chosen, drop(data$support))

# Plot
# ---------------------------------------
plot(ctlars_obj$get_ssr(), type = "l", xlab = "Lars Step", ylab = "SSR")
grid()

plot(ctlars_obj$get_r2(), type = "l", xlab = "Lars Step", ylab = "R2")
grid()
