# R/generate_ccd_data
#
#' @title
#' Generate Complex Covariate-Response Data
#'
#' @description
#' This function generates a dataset with complex covariates and a complex response
#' for use in simulations or modeling.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of predictors.
#' @param mean_xr Numeric. Mean of the real part of the predictors. Default is 0.
#' @param sd_xr Numeric. Standard deviation of the real part of the predictors. Default is 1.
#' @param mean_xi Numeric. Mean of the imaginary part of the predictors. Default is 0.
#' @param sd_xi Numeric. Standard deviation of the imaginary part of the predictors. Default is 1.
#' @param beta_cardinality Integer. Number of non-zero elements in the true beta vector. Default is 5.
#' @param noise_power Numeric. The power of the noise. Default is 1. Overridden if set_snr is TRUE.
#' @param set_snr Logical. If TRUE, the signal-to-noise ratio (SNR) will be set according to the specified value. Default is TRUE.
#' @param snr_is_linear Logical. If TRUE, the SNR value is interpreted as a linear value. If FALSE, the SNR value is interpreted in decibels. Default is TRUE.
#' @param snr_val_linear Numeric. The desired signal-to-noise ratio (SNR) as a linear value. Default is 2.
#' @param snr_val_db Numeric. The desired signal-to-noise ratio (SNR) in decibels. Default is 10.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{X}{Matrix. The generated complex predictor matrix.}
#'   \item{y}{Complex vector. The generated response vector.}
#'   \item{support}{Integer vector. Indices of the non-zero elements in the true beta vector.}
#'   \item{snr_linear}{Numeric. The signal-to-noise ratio used in the data generation process.}
#' }
#'
#' @importFrom stats rnorm
#'
#' @examples
#' # Generate a dataset with 100 observations and 10 predictors
#' data <- generate_ccg_data(n = 100, p = 10)
#' str(data)
#'
#' @export
# -----------------------------------------------------------------------------
generate_ccg_data <- function(n,
                              p,
                              mean_xr = 0,
                              sd_xr = 1,
                              mean_xi = 0,
                              sd_xi = 1,
                              beta_cardinality = 5,
                              noise_power = 1,
                              set_snr = TRUE,
                              snr_is_linear = TRUE,
                              snr_val_linear = 2,
                              snr_val_db = 10) {

  # Generate predictor matrix
  x <- matrix(data = complex(
    real = stats::rnorm(n * p, mean = mean_xr, sd = sd_xr),
    imaginary = stats::rnorm(n * p, mean = mean_xi, sd = sd_xi),
  ), nrow = n, ncol = p)

  # Define the true relationship
  beta_true_support <- sample(p, beta_cardinality)
  beta <- complex(real = rep(0, times = p),
                  imaginary = rep(0, times = p))
  beta[beta_true_support] <- complex(real = 1, imaginary = 1)

  # Calculate true response (without noise)
  y_true <- x[, beta_true_support] %*% beta[beta_true_support]

  if (set_snr) {
    if (!snr_is_linear) {
      snr_val_linear <- 10^(snr_val_db / 10)
    }
    # Calculate signal power
    signal_power <- mean(Mod(y_true) ** 2)

    # Calculate noise power
    noise_power <- signal_power / snr_val_linear
  }

  # Add complex circularly centered white noise
  z <- complex(real = stats::rnorm(n, mean = 0, sd = sqrt(noise_power / 2)),
               imaginary = stats::rnorm(n, mean = 0, sd = sqrt(noise_power / 2)))

  # Combine noise with signal to get observed response
  y <- y_true + z

  # Output
  data <- list()
  data$X <- x
  data$y <- y
  data$support <- beta_true_support
  data$snr_linear <- snr_val_linear

  return(data)
}
