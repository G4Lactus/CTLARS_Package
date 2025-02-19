rm(list = ls())

generate_simdata <- TRUE
if (generate_simdata) {
  n <- 100
  p <- 150
  beta_cardinality <- 3
  #set.seed(24)

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
                            snr_val_linear = 1,
                            snr_val_db = 10)

} else {
  library(R.matlab)
  import_data_lst <- c(
    "WGNData_rnd2.mat", # n = 10, p = 15
    "WGNData_rnd22.mat" # n = 100, p = 150
  )
  data <- R.matlab::readMat(import_data_lst[2])
}


# Append dummies
# ----------------------------------------
#set.seed(789)
update_data <- function(data_lst, num_dummies = p) {
  data_lst$p <- ncol(data_lst$X)
  data_lst$n <- nrow(data_lst$X)

  data_lst$num_dummies <- num_dummies
  data_lst$X <- cbind(data_lst$X,
                      matrix(complex(real = stats::rnorm(data_lst$n * num_dummies, sd = 1),
                                     imag = stats::rnorm(data_lst$n * num_dummies, sd = 1)),
                             # circle border
                             #exp(1i * stats::rnorm(data_lst$n*num_dummies, 0, 2 * pi)),
                             nrow = data_lst$n, ncol = num_dummies)
  )

  return(data_lst)
}
data <- update_data(data, num_dummies = ncol(data$X))


# Run complex terminating lars algorithm
# ----------------------------------------
cat("Start Complex Terminating Lars\n")
tic_init <- tictoc::tic()
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies,
                         verbose = TRUE)
toc_init <- tictoc::toc()
tic_run <- tictoc::tic()
ctlars_obj$execute_clars_step(t_stop = 1, early_stop = TRUE, use_chol = TRUE)
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

lars_path <- ctlars_obj$get_beta_history()
eps <- .Machine$double.eps
dummy_num_path <- colSums(
  matrix(
    abs(lars_path[seq.int(from = p+1, to=p+data$num_dummies), ]) > .Machine$double.eps,
    nrow = data$num_dummies,
    ncol = ncol(lars_path)
  )
  )
T_stop <- 1
K <- 1
phi_T_mat <- matrix(0, nrow = p, ncol = T_stop)
for (c in seq(T_stop)) {
  if (!any(dummy_num_path == c)) {
    ind_sol_path <- length(dummy_num_path)
    warning(
      paste(
        "For T_stop = ",
        c,
        " LARS is running until k = min(n, p) and stops there before selecting ",
        c,
        " dummies.",
        sep = ""
      )
    )
  } else {
    ind_sol_path <- which(as.numeric(dummy_num_path) == c)[1]
  }
  phi_T_mat[, c] <- (1 / K) * (abs(lars_path[1:p, ind_sol_path]) > eps)
}



# Warm Start
# -------------
ctlars_obj$execute_clars_step(t_stop = 5)
cat("Dummies in proposal: ", ctlars_obj$get_active_dummies(), "\n")
plot(ctlars_obj)
ctlars_obj$get_active_set()

## Try a second example with t_stop = 5
## ---------------------------------------
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies)
ctlars_obj$execute_clars_step(t_stop = 5, early_stop = TRUE)

cat("Active index set: ", ctlars_obj$get_active_set(), "\n")
plot(ctlars_obj)


# Full Lars path
# ----------------------------------------
run_full_path <- FALSE
if (run_full_path) {
  ctlars_obj <- ctlars$new(data$X,
                           data$y,
                           has_intercept = FALSE,
                           num_dummies = data$num_dummies)

  ctlars_obj$execute_clars_step(t_stop = 1, early_stop = FALSE)
  plot(ctlars_obj)

  ## Evaluate FDP and TPP
  ## ----------------------------------------
  chosen <- ctlars_obj$get_active_set()
  ## remove dummies
  chosen <- setdiff(
    chosen, ctlars_obj$get_active_dummies()
  )
  cat("FDP: ", FDP(chosen, drop(data$support)), "\n")
  cat("TPP: ", TPP(chosen, drop(data$support)), "\n")
  #
  ## Plot
  ## ---------------------------------------
  plot(ctlars_obj$get_ssr(), type = "l", xlab = "Lars Step", ylab = "SSR")
  grid()

  plot(ctlars_obj$get_r2(), type = "l", xlab = "Lars Step", ylab = "R2")
  grid()
}
