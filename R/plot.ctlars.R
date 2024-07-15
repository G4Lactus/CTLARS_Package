#' @title
#' Plots the CT-LARS solution path.
#'
#' @description
#' Plots the CT-LARS solution path.
#' (see [ctlars] for details.)
#'
#' @param x Object of the class ctlars See [ctlars] for details.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-axis.
#' @param include_dummies Logical. If TRUE solution paths of dummies are added to the plot.
#' @param actions Logical. If TRUE axis above plot with indices of added variables
#' (Dummies represented by 'D') along the solution path is added.
#' @param col_selected Color of lines corresponding to selected variables.
#' @param col_dummies Color of lines corresponding to included dummies.
#' @param lty_selected Line type of lines corresponding to selected variables.
#' See [par] for more details.
#' @param lty_dummies Line type of lines corresponding to included dummies.
#' See [par] for more details.
#' @param legend_pos Legend position. See [xy.coords] for more details.
#' @param ... Ignored. Only added to keep structure of generic [plot] function.
#'
#' @return Plots the ctlars solution path stored in a ctlars object.
#'
#' @importFrom stats rnorm
#' @importFrom graphics matplot axis abline mtext legend
#' @import methods
#'
#' @export
#'
#' @seealso [ctlars], [plot], [par], and [xy.coords].
#'
#' @examples
#' set.seed(42)
# -----------------------------------------------------------------------------
plot.ctlars <- function(x,
                        xlab = "# Included dummies",
                        ylab = "Mod(Coefficients)",
                        include_dummies = TRUE,
                        actions = TRUE,
                        col_selected = "black",
                        col_dummies = "red",
                        lty_selected = "solid",
                        lty_dummies = "dashed",
                        legend_pos = "topleft",
                        ...) {

  if (!methods::is(object = x, class2 = "ctlars")) {
    stop("'x' must be of class ctlars.")
  }

  T_stop <- length(x$get_active_dummies())
  num_dummies <- x$get_num_dummies()
  var_select_path <- x$get_active_set()

  beta_path <- x$get_beta_history()
  beta_path <- cbind(rep(0, times = nrow(beta_path)), beta_path)
  colnames(beta_path) <- NULL
  beta_path <- t(beta_path)
  beta_path_mod <- Mod(beta_path)
  p <- x$get_p()
  dummies_path <- which(var_select_path > p) + 1
  dummies_path_labels <- seq(T_stop)

  graphics::matplot(
    beta_path_mod[, seq(1, p)],
    col = col_selected,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    lty = lty_selected,
    xaxt = "n"
  )

  graphics::axis(
    side = 1,
    at = dummies_path,
    labels = dummies_path_labels,
    ...
  )

  graphics::abline(
    v = dummies_path,
    col = col_dummies,
    lty = 1,
    lwd = 1.3
  )

  # Add dummies solution path to plot
  if (include_dummies) {
    graphics::matlines(
      beta_path_mod[, seq(p + 1, p + num_dummies)],
      col = "red",
      type = "l",
      lty = "dashed"
    )
  }

  # Add axis above plot to indicate index of added or removed variables
  # (added dummies are indicated with 'D')
  if (actions) {
    var_select_path_positions <- seq(2, length(var_select_path) + 1)
    var_select_path_labels <- var_select_path
    var_select_path_labels[var_select_path_labels > p] <- "D"

    graphics::axis(
      side = 3,
      at = var_select_path_positions,
      labels = var_select_path_labels
    )

    graphics::mtext(
      "Index of selected variables (D indicates an included dummy)",
      side = 3,
      line = 3
    )

    graphics::abline(
      v = var_select_path_positions,
      col = "gray",
      lty = 6
    )
  }

  # Add legend to plot if active variables and dummies are plotted
  if (include_dummies && !is.null(legend_pos)) {
    graphics::legend(
      "topleft",
      legend = c("Original variables", "Dummies"),
      col = c("black", "red"),
      lty = c("solid", "dashed"),
      lwd = rep(1, times = 2)
    )
  }



}
