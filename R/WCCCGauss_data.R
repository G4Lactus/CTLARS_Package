# R/WCCCGauss_data
#
#' WCCCGauss Data
#'
#' A dataset containing complex numbers and related parameters for an example
#' of the CTLARS algorithm.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{X}{A complex matrix with n = 100 rows and p = 150 columns, containing complex numbers.}
#'   \item{beta}{A complex vector of length 150.}
#'   \item{n}{An integer representing the number of observations, which is 100.}
#'   \item{p}{An integer representing the number of variables, which is 150.}
#'   \item{support}{An integer vector indicating the support indices.}
#'   \item{y}{A complex matrix with one column and 100 rows.}
#'   \item{beta_card}{An integer representing the cardinality of beta, which is 5.}
#'   \item{rnd_state}{An integer representing the random state, which is 22.}
#' }
#'
#' @usage data(WCCCGauss_data)
#' @docType data
#'
"WCCCGauss_data"
