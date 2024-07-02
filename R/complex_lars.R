complex_lars <- function(x,
                         y,
                         desired_model_size = Inf,
                         has_intercept = TRUE,
                         standardize = TRUE,
                         verbose = TRUE,
                         tol = .Machine$double.eps) {

  # measure input values
  n_rows <- nrow(x)  # Number of rows in x
  p_cols <- ncol(x)  # Number of columns in x
  len_y <- length(y)

  # Input fidelity checks
  check_data_fidelity(len_y, n_rows, desired_model_size)

  # effective n reduced by 1 if intercept is TRUE
  effective_n <- set_effective_n(n_rows, has_intercept)

  # max Lars steps rule of thumb (according to Hastie et. al.)
  max_steps <- set_max_lars_steps(p_cols, effective_n)

  # Initialize Lars step counter
  lars_step_k <- 1

  # Data pre-processing
  # ----------------------
  normalize_data(cenv = environment())

  # Container initialization
  # --------------------------
  # Initialize the vector of regression coefficients
  beta_zm <- numeric(p_cols)

  # Initialize the tracking matrix for predicted values of y
  y_hat <- matrix(0, nrow = n_rows, ncol = 1)

  # Initialize the vector of regression coefficient updates
  gamma_hat <- numeric(p_cols)

  # Initialize the vector of active covariate indices
  index_actives <- c()

  # Set up list to track the history of indices
  coef_index_track <- list()

  # Integer representing the current active number of covariates
  p_active <- 0

  # Interger representing the current inactive number of covariates
  p_inactive <- p_cols

  # sequence of available predictors
  p_possible <- seq(p_cols)

  # intialize residuals
  res <- y

  # initialize summed squared response and summed squared residuals
  ssy <- sum(Mod(y)**2)

  # Initialize statistical measures
  ssr <- ssy
  r2 <- c(0)

  # Initialize new index
  new_index <- NULL

  # Modified tolerance variable
  tolc <- 1 - tol

  # conjugate transpose of x
  x_ct <- Conj(t(x))

  while (p_active < effective_n &&
           p_inactive > 0 &&
           lars_step_k < max_steps &&
           length(new_index) + p_active <= desired_model_size) {

    # Step (2.1)
    # Find correlation of each covariate with the current residual
    c_hat <- x_ct %*% res

    # Step (2.2) Find maximum absolute current correlation
    c_tild <- max(Mod(c_hat))

    # 1st iteration: select the indices of the maximum absolute current
    # correlations to the active set
    if (lars_step_k == 1) {
      # Find index of every correlation that is within tolerance
      # of the max correlation
      c_max_ind <- which(Mod(c_hat) > tolc * c_tild)
      new_index <- c_max_ind
    }
    progress_report(verbose, lars_step_k)

    # Add the appropriate indices to the active set
    index_actives <- c(index_actives, new_index)
    coef_index_track[[lars_step_k]] <- index_actives

    # Calculate the number of indices in the active set
    p_active <- length(index_actives)

    # Update the set of inactive indices
    compl_set <- setdiff(p_possible, index_actives)

    # Calculate the number of remaining inactive candidate covariates
    p_inactive <- p_cols - p_active

    # Create a vector containing only the active correlations
    # (each entry should be equal in magnitude)
    c_hat_s <- c_hat[index_actives]

    # Calculate the complex signum values of each active correlation
    sign_val <- complex_signum(c_hat_s)

    # Initialize matrix to store active, signum-aligned covariates
    sigaligned_x_s <- matrix(0, nrow = n_rows, ncol = p_active)
    # signum align the matrix active set candidates
    for (i in seq_len(p_active)) {
      sigaligned_x_s[, i] <- x[, index_actives[i]] * sign_val[i]
    }

    # Step (2.4)
    # Create Gramian matrix
    g_mat <- Conj(t(sigaligned_x_s)) %*% sigaligned_x_s

    # Invert Gramian matrix
    inv_g_mat <- compute_gram_inverse(g_mat, lars_step_k, tol)

    # one vector
    one <- rep(1, times = p_active)
    # norm factor for equi-angle condition
    invg_one <- inv_g_mat %*% one
    l_a <- Re((t(one) %*% invg_one)**(-0.5))
    # weights for equi-angle direction
    w_a <- c(l_a) * invg_one
    # equi-angle direction
    u_a <- sigaligned_x_s %*% w_a

    # step (2.5)
    # direction of travel for all covariates
    g <- x_ct %*% u_a

    # Step (2.6) find the length gamma_hat to travel u_a
    gamma_hat <- update_gamma_hat(p_active, p_cols, lars_step_k,
                                  c_tild, l_a, p_inactive, compl_set,
                                  tol, c_hat, g, gamma_hat,
                                  cenv = environment())

    # update coefficients
    beta_zm <- update_beta_zm(p_cols, gamma_hat,
                              lars_step_k, sign_val,
                              w_a, index_actives,
                              cenv = environment())

    # update predictions
    y_hat <- y_hat + gamma_hat[lars_step_k] * u_a

    # increment clars step
    lars_step_k <- lars_step_k + 1

    # update residual
    res <- y - y_hat

    # Update statistics
    update_statistics(cenv = environment())

    # End of CLARS step
  }

  output <- list()
  output$index_actives <- index_actives
  output$coef_index_track <- coef_index_track
  output$y_hat <- y_hat
  output$beta_zm <- Conj(t(beta_zm))
  output$ssr <- ssr
  output$r2 <- r2

  return(output)
}


scale_x <- function(x, num_rows, num_cols) {

  # center values
  x <- scale(x, center = TRUE, scale = FALSE)

  # determine column standard deviations
  x_cstds <- apply(x, 2, function(xs) {
    return(sqrt(sum(Mod(xs) ** 2)  / (num_rows - 1)))
  })

  # normalize
  for (colX in seq_len(num_cols)) {
    x[, colX] <- x[, colX] / x_cstds[colX]
  }
  return(x)
}


set_effective_n <- function(n_rows, has_intercept) {
  effective_n <- n_rows
  if (has_intercept) {
    effective_n <- effective_n - 1
  }
  return(effective_n)
}


check_data_fidelity <- function(len_y, n_rows, desired_model_size) {
  if (len_y != n_rows) {
    stop("Covariates must be the same length as the data.")
  }
  if (!is.numeric(desired_model_size) || desired_model_size <= 0) {
    stop("Your desired model size must be numeric and at least 1.")
  }
}


set_max_lars_steps <- function(p_cols, effective_n) {
  if (p_cols < effective_n) {
    return(8 * p_cols)
  }
  return(8 * effective_n)
}


progress_report <- function(verbose, lars_step_k) {
  if (verbose) {
    cat("complex LARS iteration", lars_step_k, "\n")
  }
}



update_beta_zm <- function(p_cols, gamma_hat,
                           lars_step_k, sign_val,
                           w_a, index_actives,
                           cenv) {
  beta_tmp <- numeric(p_cols)
  if (lars_step_k != 1) {
    # note: t() here omitted due to dimensionless beta_tmp
    beta_tmp[index_actives] <-
      Conj(cenv$beta_zm[lars_step_k - 1, index_actives]) +
      gamma_hat[lars_step_k] * (sign_val * w_a)
    return(rbind(cenv$beta_zm, Conj(t(beta_tmp))))
  } else {
    # beta_tmp[index_actives] <- Conj(t(beta_zm[index_actives])) +
    beta_tmp[index_actives] <- gamma_hat[lars_step_k] * (sign_val * w_a)
    return(Conj(t(beta_tmp)))
  }
}


update_gamma_hat <- function(p_active, p_cols, lars_step_k,
                             c_tild, l_a, p_inactive, compl_set,
                             tol, c_hat, g, gamma_hat, cenv) {

  if (p_active == p_cols) {  # (compl_set is empty)
    # move to the least squares projection
    gamma_hat[lars_step_k] <- c_tild / l_a

  } else { # (compl_set is not empty)
    gamma_candids <- matrix(0, nrow = p_inactive, ncol = 2)
    for (i in seq_len(p_inactive)) {
      ii <- compl_set[i]
      a <- Mod(g[ii]) ** 2 - (l_a ** 2) + tol
      b <- -2 * ((Re(Conj(t(g[ii])) %*% c_hat[ii])) - (c_tild * l_a))
      c <- (Mod(c_hat[ii]) ** 2) - (c_tild ** 2)
      d <- (b ** 2) - (4 * a * c) + tol
      if (d < tol) {
        d <- 1e-2 * tol
      }
      gamma_candids[i, 1] <- (-b + sqrt(d)) / (2 * a)
      gamma_candids[i, 2] <- (-b - sqrt(d)) / (2 * a)
    }
    gamma_candids[gamma_candids <= 0] <- Inf
    min_gam_val_ind <- apply(gamma_candids, 1, function(row_x) {
      return(c(min(row_x), which.min(row_x)))
    })
    min_gam_1 <- min_gam_val_ind[1, ]
    min_ind_1 <- min_gam_val_ind[2, ]
    min_gam_2 <- min(min_gam_1)
    min_ind_2 <- which.min(min_gam_1)
    min_ind_1 <- min_ind_1[min_ind_2]
    min_ind_2 <- which(Mod(min_gam_1 - min_gam_2) < tol * min_gam_2)

    cenv$new_index <- compl_set[min_ind_2]  # global update
    gamma_hat[lars_step_k] <- min_gam_2
  }
  return(gamma_hat)
}


complex_signum <- function(z) {

  return(
    sapply(z, function(z_val) {
      if (z_val == 0) {
        return(0)
      } else {
        return(z_val / Mod(z_val))
      }
    })
  )
}


update_statistics <- function(cenv) {
  vssr <- sum(Mod(cenv$res)**2)
  cenv$ssr <- c(cenv$ssr, vssr)
  cenv$r2 <- c(cenv$r2, 1 - vssr / cenv$ssy)
}


compute_gram_inverse <- function(gram_mat, lars_step_k, tol_val) {
  tryCatch({
    # Attempt to compute the inverse with the default tolerance level
    inv_g_mat <- Matrix::solve(gram_mat, tol = 1e-2 * tol_val)
    return(inv_g_mat)
  }, error = function(e) {
    # Check if the error message indicates computational singularity
    if (grepl("system is computationally singular", conditionMessage(e))) {
      # If singularity is detected, reduce the tolerance level and retry
      reduced_tol <- tol_val ** lars_step_k
      inv_g_mat <- Matrix::solve(gram_mat, tol = reduced_tol)
      return(inv_g_mat)
    } else {
      # If it's a different error, re-throw the error
      stop(e)
    }
  })
}


normalize_data <- function(cenv) {
  if (cenv$standardize && cenv$n_rows > 1) {
    # standardize matrix X
    cenv$x <- scale_x(cenv$x, num_rows = cenv$n_rows, num_cols = cenv$p_cols)

    # centralize vector y
    cenv$y <- scale(cenv$y, center = TRUE, scale = FALSE)
  }
}
# ---------------------------------------