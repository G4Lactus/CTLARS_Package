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
                          snr_val_linear = 0.001
                          )

l <- 5
data$X_augmented <- cbind(data$X,
                          matrix(complex(
                            real = stats::rnorm(n * l * p),
                            imaginary = stats::rnorm(n * l * p)),
                            nrow = n, ncol = l * p))


ctlars_obj <- ctlars::ctlars$new(data$X_augmented, data$y,
                   has_intercept = FALSE,
                   standardize = TRUE,
                   num_dummies = p, verbose = TRUE)

ctlars_obj$execute_clars_step(4)
cat("True support: ", data$support, "\n")
cat("Active set: ", ctlars_obj$get_active_set(), "\n")

