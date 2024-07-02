clars <- R6::R6Class(
  "clars",
  public = list(
    # Public attributes
    # --------------------
    x = NULL,
    y = NULL,
    has_intercept = NULL,
    standardize = NULL,
    verbose = NULL,

    # Public methods
    # ---------------------
    # constructor method
    initialize = function(
      x,
      y,
      has_intercept = TRUE,
      standardize = TRUE,
      verbose = TRUE,
      tol = .Machine$double.eps
    ) {
      # set attributes
      self$x <- x
      self$y <- y
      self$has_intercept <- has_intercept
      self$standardize <- standardize
      self$verbose <- verbose
      private$tol <- tol

      # prepare inputs
      private$initialize_values()
    },

    # execute clars step
    execute_clars_step = function(desired_model_size = Inf) {
      private$check_model_size(desired_model_size)

      while (private$p_active < private$effective_n &&
          private$p_inactive > 0 &&
          private$lars_step_k < private$max_steps &&
          length(private$new_index) + private$p_active <= desired_model_size
      ) {

        # Step (2.1)
        # Find correlation of each covariate with the current residual
        c_hat <- private$x_ct %*% private$res

        # Step (2.2)
        # Find maximum absolute current correlation
        c_tild <- max(Mod(c_hat))

        # Step (2.3)
        # 1st iteration: select the indices of the maximum absolute current
        # correlations to the active set
        if (private$lars_step_k == 1) {
          # Find index of every correlation that is within tolerance
          # of the max correlation
          c_max_ind <- which(abs(c_hat) > private$tolc * c_tild)
          private$new_index <- c_max_ind
        }
        private$progress_report()

        # Add the appropriate indices to the active set
        private$index_actives <- c(private$index_actives, private$new_index)
        private$coef_index_track[[private$lars_step_k]] <- private$index_actives

        # Calculate the number of indices in the active set
        private$p_active <- length(private$index_actives)

        # Update the set of inactive indices
        compl_set <- setdiff(private$p_possible, private$index_actives)

        # Calculate the number of remaining inactive candidate covariates
        private$p_inactive <- private$p_cols - private$p_active

        # Create a vector containing only the active correlations
        # (each entry should be equal in magnitude)
        c_hat_s <- c_hat[private$index_actives]

        # Calculate the complex signum values of each active correlation
        sign_val <- private$complex_signum(c_hat_s)

        # Initialize matrix to store active, signum-aligned covariates
        sigaligned_x_s <- matrix(0,
                                 nrow = private$n_rows,
                                 ncol = private$p_active)

        # signum align the matrix active set candidates
        for (i in seq_len(private$p_active)) {
          sigaligned_x_s[, i] <- self$x[, private$index_actives[i]] *
            sign_val[i]
        }

        # Step (2.4)
        # Create Gramian matrix
        g_mat <- Conj(t(sigaligned_x_s)) %*% sigaligned_x_s

        # Invert Gramian matrix
        inv_g_mat <- private$compute_gram_inverse(g_mat)

        # one vector
        one <- rep(1, times = private$p_active)
        # norm factor for equi-angle condition
        invg_one <- inv_g_mat %*% one
        l_a <- Re((t(one) %*% invg_one)**(-0.5))
        # weights for equi-angle direction
        w_a <- c(l_a) * invg_one
        # equi-angle direction
        u_a <- sigaligned_x_s %*% w_a

        # step (2.5)
        # direction of travel for all covariates
        g <- private$x_ct %*% u_a

        # Step (2.6) find the length gamma_hat to travel u_a
        private$update_gamma_hat(
          c_tild,
          l_a,
          compl_set,
          c_hat,
          g
        )

        # update coefficients
        private$update_beta_zm(sign_val, w_a)

        # update predictions
        private$y_hat <- private$y_hat +
          private$gamma_hat[private$lars_step_k] * u_a

        # increment clars step
        private$lars_step_k <- private$lars_step_k + 1

        # update residual
        private$res <- self$y - private$y_hat

        # update statistics
        private$update_statistics()

        # End of CLARS step
      }
    },

    get_active_set = function() {
      return(private$index_actives)
    },

    get_beta_hat = function() {
      return(Conj(t(private$beta_zm))[, private$lars_step_k])
    },

    get_beta_history = function() {
      return(Conj(t(private$beta_zm)))
    },

    get_y_hat = function() {
      return(private$y_hat)
    },

    get_ssr = function() {
      return(private$ssr)
    },

    get_r2 = function() {
      return(private$r2)
    }
  ),
  private = list(
    # private attributes
    # -------------------
    n_rows = NULL,
    p_cols = NULL,
    len_y  = NULL,
    new_index = NULL,
    tol = NULL,
    tolc = NULL,
    effective_n = NULL,
    max_steps = NULL,
    beta_zm = NULL,
    res = NULL,
    ssy = NULL,
    ssr = NULL,
    r2 = NULL,
    y_hat = NULL,
    gamma_hat = NULL,
    index_actives = NULL,
    coef_index_track = NULL,
    p_active = NULL,
    p_inactive = NULL,
    p_possible = NULL,
    lars_step_k = NULL,
    x_ct = NULL,

    # private methods
    # ---------------
    # initialize method
    initialize_values = function() {
      private$n_rows <- nrow(self$x)
      private$p_cols <- ncol(self$x)
      private$len_y <- length(self$y)

      # Input data fidelity checks
      private$check_data_fidelity()

      # effective n reduced by 1 if intercept is TRUE
      private$set_effective_n()

      # max CLars steps rule of thumb (according to Hastie et. al.)
      private$set_max_clars_steps()

      # initialize CLars step counter
      private$lars_step_k <- 1

      # Data preprocessing
      # -------------------
      private$normalize_data()

      # Container initialization
      # ------------------------
      # Initialize the vector of regression coefficients
      private$beta_zm <- numeric(private$p_cols)

      # Initialize the tracking matrix for predicted values of y
      private$y_hat <- matrix(data = 0, nrow = private$n_rows, ncol = 1)

      # Initialize the vector of active covariate indices
      private$gamma_hat <- numeric(private$p_cols)

      # Initialize the vector of active covariate indices
      private$index_actives <- c()

      # Set up a list to track the history of indices
      private$coef_index_track <- list()

      # Integer representing the current active number of covariates
      private$p_active <- 0

      # Integer representing the current inactive number of covariates
      private$p_inactive <- private$p_cols

      # Sequence of available predictors
      private$p_possible <- seq(private$p_cols)

      # initialize residuals
      private$res <- self$y

      # initialize summed squared response and summed squared residuals
      private$ssy <- sum(Mod(self$y)**2)

      # Initialize statistical measures
      private$ssr <- private$ssy
      private$r2 <- c(0)

      # Modified tolerance variable
      private$tolc <- 1 - private$tol

      # store conjugate transpose of x
      # can also be placed in execute_clars_step, then temporary
      private$x_ct <- Conj(t(self$x))
    },

    check_data_fidelity = function() {
      if (private$len_y != private$n_rows) {
        stop("The number of observations `y` must match `n`.")
      }
    },

    check_model_size = function(desired_model_size) {
      if (!is.numeric(desired_model_size) || desired_model_size <= 0) {
        stop("Your `desired_model_size` must be numeric and at least 1.")
      }
    },

    set_effective_n = function() {
      n <- private$n_rows
      if (self$has_intercept && n > 1) {
        private$effective_n <- n - 1
      } else {
        private$effective_n <- n
      }
    },

    set_max_clars_steps = function() {
      p_cols <- private$p_cols
      n_eff <- private$effective_n
      if (p_cols < n_eff) {
        private$max_steps <- 8 * p_cols
      } else {
        private$max_steps <- 8 * n_eff
      }
    },

    normalize_data = function() {
      if (self$standardize && private$n_rows > 1) {
        # standardize matrix X
        private$scale_x()
        # centralize vector y
        self$y <- scale(self$y, center = TRUE, scale = FALSE)
      }
    },

    scale_x = function() {
      x <- self$x
      # center values
      x <- scale(x, center = TRUE, scale = FALSE)
      n <- private$n_rows
      if (n < 2) {
        denum <- n
      } else {
        denum <- n - 1
      }

      # normalize
      for (colX in seq_len(private$p_cols)) {
        col <- x[, colX]
        col_std <- sqrt(sum(Mod(col) ** 2) / denum)
        x[, colX] <- col / col_std
      }
      self$x <- x
    },

    progress_report = function() {
      if (self$verbose) {
        cat("complex LARS iteration", private$lars_step_k, "\n")
      }
    },

    complex_signum = function(z) {
      return(
        ifelse(z == 0, 0, z / Mod(z))
      )
    },

    compute_gram_inverse = function(gram_mat) {
      tryCatch({
        # Attempt to compute the inverse with the default tolerance level
        return(Matrix::solve(gram_mat, tol = 1e-2 * private$tol_val))
      }, error = function(e) {
        # Check if the error message indicates computational singularity
        if (grepl("System is computationally singular.", conditionMessage(e))) {
          # If singularity is detected, reduce the tolerance level and retry
          reduced_tol <- private$tol_val ** private$lars_step_k
          return(Matrix::solve(gram_mat, tol = reduced_tol))
        } else {
          # If it's a different error, re-throw the error
          stop(e)
        }
      })
    },

    update_beta_zm = function(
      sign_val,
      w_a
    ) {
      beta_tmp <- numeric(private$p_cols)
      beta_zm <- private$beta_zm
      if (private$lars_step_k != 1) {
        # note: t() here omitted due to dimensionless beta_tmp
        beta_tmp[private$index_actives] <-
          Conj(beta_zm[private$lars_step_k - 1,
                       private$index_actives]) +
          private$gamma_hat[private$lars_step_k] * (sign_val * w_a)
      } else {
        beta_tmp[private$index_actives] <-
          private$gamma_hat[private$lars_step_k] * (sign_val * w_a)
      }
      private$beta_zm <- rbind(beta_zm, Conj(t(beta_tmp)))
    },

    update_gamma_hat = function(
      c_tild,
      l_a,
      compl_set,
      c_hat,
      g
    ) {
      tol <- private$tol
      if (private$p_active == private$p_cols) {  # (compl_set is empty)
        # move to the least squares projection
        private$gamma_hat[private$lars_step_k] <- c_tild / l_a

      } else { # (compl_set is not empty)
        gamma_candids <- matrix(0, nrow = private$p_inactive, ncol = 2)
        mod_g_sq <- Mod(g) ** 2
        mod_c_sq <- Mod(c_hat) ** 2
        for (i in seq_len(private$p_inactive)) {
          ii <- compl_set[i]
          a <- mod_g_sq[ii] - (l_a ** 2) + tol
          b <- -2 * ((Re(Conj(t(g[ii])) %*% c_hat[ii])) - (c_tild * l_a))
          c <- mod_c_sq[ii] - (c_tild ** 2)
          d <- max((b ** 2) - (4 * a * c) + tol, 1e-2 * tol)
          gamma_candids[i, 1] <- (-b + sqrt(d)) / (2 * a)
          gamma_candids[i, 2] <- (-b - sqrt(d)) / (2 * a)
        }
        gamma_candids[gamma_candids <= 0] <- Inf
        min_gam_1 <- apply(gamma_candids, 1, min)
        min_gam_2 <- min(min_gam_1)
        min_ind_2 <- which(Mod(min_gam_1 - min_gam_2) < tol * min_gam_2)

        # ------------
        private$new_index <- compl_set[min_ind_2]
        private$gamma_hat[private$lars_step_k] <- min_gam_2
      }
    },

    update_statistics = function() {
      vssr <- sum(Mod(private$res)**2)
      private$ssr <- c(private$ssr, vssr)
      private$r2 <- c(private$r2, 1 - vssr / private$ssy)
    }
  )
  # ---------------
  # End of class
)