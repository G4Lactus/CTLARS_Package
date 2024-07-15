# R/WCCCGauss_data_small
#
#' WCCCGauss Data small
#'
#' A dataset containing complex numbers and related parameters for testing
#' the CTLARS algorithm.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{X}{A complex matrix with 10 rows and 15 columns, containing complex numbers.}
#'   \item{beta}{A complex vector of length 15.}
#'   \item{n}{An integer representing the number of observations, which is 10.}
#'   \item{p}{An integer representing the number of variables, which is 15.}
#'   \item{support}{An integer vector indicating the support indices.}
#'   \item{y}{A complex matrix with one column and 10 rows.}
#'   \item{beta_card}{An integer representing the cardinality of beta, which is 5.}
#'   \item{rnd_state}{An integer representing the random state, which is 2.}
#' }
#'
#' @usage data(WCCCGauss_data_small)
#' @docType data
#'
"WCCCGauss_data_small"
