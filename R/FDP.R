# R/FDP
#
#' @title
#' False Detection Proportion
#'
#' @description
#' The function returns the empirical false detection proportion based on a
#' vector of selected variable indices and a vector of reference indices.
#'
#' @name
#' FDP
#'
#' @usage
#' FDP(selected_vars, true_actives)
#'
#' @param selected_vars a vector with estimated support indices.
#' @param true_actives a vector with true support indices.
#'
#' @return
#' Returns the false detection proportion as a number.
#'
#' @export
#'
#' @examples
#' set.seed(42)
# -----------------------------------------------------------------------------
FDP <- function(selected_vars, true_actives) {

  # Error control
  if (!is.vector(selected_vars)) {
    stop("'selected_vars' must be a vector.")
  }

  if (!is.numeric(selected_vars)) {
    stop("'selected_vars' only allows numerical values.")
  }

  if (anyNA(selected_vars)) {
    stop("'selected_vars' contains NAs. Please remove or impute them before proceeding.")
  }

  if (!is.vector(true_actives)) {
    stop("'true_actives' must be a vector.")
  }

  if (!is.numeric(true_actives)) {
    stop("'true_actives' only allows numerical values.")
  }

  if (anyNA(true_actives)) {
    stop("'true_actives' contains NAs. Please remove or impute them before proceeding.")
  }

  # Compute FDP
  return(length(setdiff(selected_vars, true_actives)) / max(1, length(selected_vars)))
}
# -----------------------------------------------------------------------------
