# R/simulate_model_TPR_FDR
#
#' @title
#' Simulate a ctlars model with True Positive Rate and False Detection Rate
#'
#' @description
#' From a statistical point of view, it is desirable to use a variable selection
#' method that allows for controlling the expected value of the FDP at a user
#' defined target level alpha in `[`0, 1`]`, while maximizing the number of selected
#' variables.
#' These type of methods exist and are called false discovery rate (FDR)-
#' controlling methods.
#' For example, the CT-Rex Selector is a fast and FDR-controlling variable/feature
#' selection framework for large-scale high-dimensional settings that relies on
#' the ctars method.
#'
#' Consider the following definitions of FDR and TPR:
#' Let AA be part of `{1, ..., p}` be the set of selected variables, and let
#' `A` in `{1, ..., p}` be the set of true active variables, `|AA|` be the cardinality
#' of `AA`, and define `1 v a := max{1, a}`, a in R.
#' Then, the false discovery rate (FDR) and the true positive rate (TPR) are
#' defined by:
#'
#' `FDR := E[FDP] := E[ |[AA\A]|/(1 v |AA|) ]`
#' `TPR := E[TPP] := E[ |[intersection(A, AA)]|/(1 v |A|) ]`
#'
#' This function simulates a Gaussian regression scenario with complex numbers,
#' applies the ctlars method for variable selection, and computes FDR and TPR
#' metrics across multiple termination steps. It generates noisy observations
#' and includes a specified number of active and dummy variables.
#'
#' @name
#' simulate_model_TPR_FDR
#'
#' @param T_vec Vector of terminating steps to conduct the simulation.
#'   Default is c(1, 2, 5, 10, 20, 50, 100).
#' @param MC Number of Monte Carlo experiments per termination step.
#'   Default is 50.
#' @param n Number of observations/rows for the Gaussian regression.
#'   Default is 100.
#' @param p Number of predictors/columns for the Gaussian regression.
#'   Default is 300.
#' @param num_actives Number of active variables used in Gaussian regression.
#'   Default is 10.
#' @param num_dummies Number of dummy variables added to the predictor matrix.
#'   Default is `p`.
#' @param noise_power Numeric value indicating the noise power of white
#'   centered circulary symmetric Gaussian noise. Default is 1.
#'   Can be computed by setting the parameter `set_snr` to `TRUE`.
#' @param set_snr Logical indicating if SNR (Signa-to-Noise Ratio) is controlled.
#'   Default is FALSE.
#' @param snr_is_linear Logical indicating if SNR (Signal-to-Noise Ratio) value
#'   is given in linear scale (`TRUE`) or dB scale (`FALSE`). Default is `TRUE`.
#' @param snr_val_linear SNR value in linear scale. Default is 2.
#' @param snr_val_db SNR value in dB scale. Default is 10.
#' @param verbose Logical indicating whether a simulation progress is displayed.
#' Default is `TRUE`.
#'
#' @seealso
#' \code{\link{ctlars}} for details on the ctlars algorithm.
#' \code{\link{FDP}} for details on calculating False Discovery Proportion.
#' \code{\link{TPP}} for details on calculating True Positive Proportion.
#'
#' @export
#'
#' @import ctlars
#' @importFrom stats rnorm
#'
#' @return
#' The function returns a list with the following components:
#' \itemize{
#'  \item{T_vec} Vector with realized terminating steps.
#'  \item{FDP} Matrix with realized False Discovery Proportions (FDP).
#'  \item{TPP} Matrix with realized True Positive Proportions (TPP).
#'  \item{FDR} Vector with realized False Discovery Rates (FDR).
#'  \item{TPR} Vector with realized True Positive Rates (TPR).
#' }
#'
#' @examples
#' set.seed(42)
#' result <- simulate_model_TPR_FDR()
#' print(result)
# -----------------------------------------------------------------------------
simulate_model_TPR_FDR <- function(T_vec = c(1, 2, 5, 10, 20, 50, 100),
                                   MC = 2,
                                   n = 10,
                                   p = 30,
                                   num_actives = 3,
                                   num_dummies = 30,
                                   noise_power = 1,
                                   set_snr = FALSE,
                                   snr_is_linear = TRUE,
                                   snr_val_linear = 2,
                                   snr_val_db = 10,
                                   verbose = TRUE
                                   ) {

  # error check
  if (length(T_vec) == 1) {
    warning(paste(
      "'T_vec' is supposed to be a vector",
      "to create meaningful simulation results."
    ))
  }
  if (!is.vector(T_vec)) {
    stop("Please choose 'T_vec' as vector of integer numeric values.")
  }
  if (!is.numeric(T_vec)) {
    stop("Please choose 'T_vec' as integer numeric.")
  }
  if (!is.numeric(MC) || MC <= 0) {
    stop("Please choose 'MC' as a positive integer")
  }
  if (!is.numeric(n) || n <= 0) {
    stop("Please choose 'n' as a positive integer.")
  }
  if (!is.numeric(p) || p <= 0) {
    stop("'p' must be a positive integer.")
  }
  if (!is.numeric(num_actives) || num_actives <= 0) {
    stop("Please choose 'num_actives' as positive integer.")
  }
  if (!is.numeric(num_dummies) || num_dummies < p) {
    stop("Please choose 'num_dummies' as positive integer.")
  }
  if (!is.logical(snr_is_linear)) {
    stop("Please choose 'snr_is_linear' as logical.")
  }
  if (!is.numeric(snr_val_linear) || snr_val_linear <= 0) {
    stop("Please choose 'snr_val_linear' as numeric.")
  }
  if (!is.numeric(snr_val_db)) {
    stop("Please choose 'snr_val_db' as numeric.")
  }


  # Setup beta vector
  beta_true_support <- sample(p, num_actives)
  beta <- complex(real = rep(0, times = p),
                  imaginary = rep(0, times = p))
  beta[beta_true_support] <- complex(real = 1, imaginary = 1)

  # Initialize results matrices
  FDP <- matrix(NA, nrow = MC, ncol = length(T_vec))
  TPP <- matrix(NA, nrow = MC, ncol = length(T_vec))

  # set seed for reproducibility
  set.seed(12345)

  # Perform Monte Carlo simulations
  for (t in seq_along(T_vec)) {
    for (mc in seq(MC)) {

      # Generate data
      # --------------------------
      # Generate predictor matrix
      x <- matrix(
        data = complex(
          real = matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p),
          imaginary = matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
        ),
        nrow = n,
        ncol = p
      )

      # Calculate true response (without noise)
      y_true <- x[, beta_true_support] %*% beta[beta_true_support]

      if (set_snr) {

        # Calculate noise power based on SNR
        if (!snr_is_linear) {
          snr_val_linear <- 10^(snr_val_db / 10)
        }

        # Calculate signal power
        signal_power <- mean(Mod(y_true) ^ 2)

        # Calculate noise power
        noise_power <- signal_power / snr_val_linear

      }

      # Add complex circularly centered white noise
      # Note: in the model theory, sd = sqrt(noise_power/2)
      z <- complex(real = stats::rnorm(n, mean = 0, sd = noise_power),
                   imaginary = stats::rnorm(n, mean = 0, sd = noise_power))

      # Combine noise with signal to get observed response
      y <- y_true + z

      # append dummies
      x <- cbind(x,
        matrix(data = complex(
          real = stats::rnorm(n * num_dummies),
          imaginary = stats::rnorm(n * num_dummies)
        ),
        nrow = n,
        ncol = num_dummies
        )
      )

      # Create ctlars object
      ctlars_obj <- ctlars$new(x,
                               y,
                               has_intercept = FALSE,
                               num_dummies = num_dummies,
                               verbose = verbose)

      # execute ctlars step
      ctlars_obj$execute_clars_step(t_stop = t)

      # extract results
      selected_vars <- ctlars_obj$get_active_set()
      selected_vars <- selected_vars[selected_vars < (p + 1)]

      # Results
      FDP[mc, t] <- FDP(selected_vars, beta_true_support)
      TPP[mc, t] <- TPP(selected_vars, beta_true_support)

    }
  }

  # Define output container
  output <- list()
  output$T_vec <- T_vec
  output$FDP <- FDP
  output$TPP <- TPP
  output$FDR <- colMeans(FDP)
  output$TPR <- colMeans(TPP)

  return(output)
}
# -----------------------------------------------------------------------------
