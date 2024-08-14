load("~/Documents/R/ctlars/R_demo_scripts/ULA_test.RData")

cat("Start Complex Terminating Lars\n")
ctlars_obj <- ctlars::ctlars$new(A,
                                 y,
                                 has_intercept = FALSE,
                                 num_dummies = ncol(A),
                                 verbose = TRUE)
ctlars_obj$execute_clars_step(t_stop = 30, early_stop = FALSE, use_chol = TRUE)
