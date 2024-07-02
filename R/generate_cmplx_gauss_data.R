# Scenarios:
# ----------------------------------
# - standard:
#   linear regression model with complex Gaussian data
#
# - int_conjs_in_X:
#   a pair of integer conjugate regressors x_j1 and x_j2 in covariates of X
#
# - conjs_in_X:
#   a pair of conjugate regressors x_j1 and x_j2 in covariates of X
#
# - int_conjs_in_support:
#   a pair of integer conjugate regressors x_j1 and x_j2 in support of beta_true
#
# - conjs_in_support:
#   a pair of conjugate regressors x_j1 and x_j2 in support of beta_true
#
# - copy_in_X:
#   a pair of identical regressors x_j1 and x_j2 in covariates of X
#
# - copy_in_support:
#   a pair of identical regressors x_j1 and x_j2 in support of beta_true
# ----------------------------------

generate_cmplx_gauss_data <- function(
  n_rows,
  p_cols,
  beta_cardinality = 5,
  scenario = "standard",
  beta_type = "complex", # "real", "imaginary",
  beta_val = "random" # "constant"
) {

  mat_x <- matrix(data = complex(real = stats::rnorm(n_rows * p_cols),
                                 imag = stats::rnorm(n_rows * p_cols)),
                  nrow = n_rows, ncol = p_cols)

  vec_beta <- complex(real = rep(0, times = p_cols),
                      imaginary = rep(0, times = p_cols))

  if (scenario == "standard") {

    true_support <- sample(p_cols, beta_cardinality)
    vec_beta <- generate_betas(vec_beta,
                               true_support,
                               beta_cardinality,
                               beta_type,
                               beta_val)

  } else if (scenario == "conj_regressors_outside") {

    conj_vars <- sample(p_cols, 2)
    mat_x[, conj_vars[2]] <- Conj(mat_x[, conj_vars[1]])

    true_support <- sample(setdiff(seq(p_cols), conj_vars),
                           beta_cardinality)
    vec_beta <- generate_betas(vec_beta,
                               true_support,
                               beta_cardinality,
                               beta_type,
                               beta_val)

  } else if (scenario == "conj_regressors_inside") {

    conj_vars <- sample(p_cols, 2)
    mat_x[, conj_vars[2]] <- Conj(mat_x[, conj_vars[1]])

    true_support <- sample(c(sample(setdiff(seq(p_cols), conj_vars),
                                    beta_cardinality - 2), conj_vars))
    vec_beta <- generate_betas(vec_beta,
                               true_support,
                               beta_cardinality,
                               beta_type,
                               beta_val)

  } else if (scenario == "copy_regressors_inside") {

    conj_vars <- sample(p_cols, 2)
    mat_x[, conj_vars[2]] <- mat_x[, conj_vars[1]]

    true_support <- sample(c(sample(setdiff(seq(p_cols), conj_vars),
                                    beta_cardinality - 2), conj_vars))
    vec_beta[true_support] <- complex(real = stats::rnorm(beta_cardinality),
                                      imag = stats::rnorm(beta_cardinality))
    vec_beta[conj_vars[2]] <- Conj(vec_beta[conj_vars[1]])

  } else {
    stop(cat("`scenario`: ", scenario, " currently not supported."))
  }

  # Form clean vec_y
  vec_y <- mat_x %*% vec_beta

  # Form noise vector
  vec_noise <- complex(real = stats::rnorm(n_rows),
                       imag = stats::rnorm(n_rows))

  # Form contaminated model vec_y
  vec_y <- vec_y + vec_noise

  # Generate output
  output <- list()
  output$X <- mat_x
  output$y <- vec_y
  output$support <- true_support
  output$beta <- vec_beta

  return(output)
}


generate_betas <- function(vec_beta,
                           support_index,
                           beta_cardinality,
                           beta_type,
                           beta_val) {

  if (beta_val == "constant") {
    vec_beta[support_index] <- complex(real = 1, imag = 1)

  } else if (beta_val == "random") {
    vec_beta[support_index] <- complex(
      real = stats::rnorm(beta_cardinality),
      imag = stats::rnorm(beta_cardinality)
    )
  } else {
    stop(cat("`beta_val`: ", beta_val, " currently not supported."))
  }
  if (beta_type == "complex") {
    return(vec_beta)
  } else if (beta_type == "real") {
    return(Re(vec_beta))
  } else if (beta_type == "imag")  {
    return(Im(vec_beta))
  }
}