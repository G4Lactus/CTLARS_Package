## ----include = FALSE----------------------------------------------------------
# Store user's options()
old_options <- options()

library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.retina = 2,
  out.width = "85%",
  dpi = 96
  # pngquant = "--speed=1"
)
options(width = 80)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  library(ctlars)
#  help(package = ctlars)
#  ?`ctlars-package`
#  ?ctlars
#  ?FDP
#  ?generate_ccg_data
#  ?plot.ctlars
#  ?simulate_model_TPR_FDR
#  ?TPP
#  ?WCCCGauss_data
#  ?WCCCGauss_data_small

## ----eval=FALSE---------------------------------------------------------------
#  # some future citation:
#  citation("ctlars")

## ----echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
library(ctlars)

generate_simdata <- FALSE
if (generate_simdata) {
  n <- 100
  p <- 150
  beta_cardinality <- 5
  set.seed(42)

  # generate some data
  data <- generate_ccg_data(n,
                            p,
                            mean_xr = 0,
                            sd_xr = 1,
                            mean_xi = 0,
                            sd_xi = 1,
                            beta_cardinality = beta_cardinality,
                            noise_power = 1,
                            set_snr = TRUE,
                            snr_is_linear = TRUE,
                            snr_val_linear = 2,
                            snr_val_db = 10)

} else {
  data("WCCCGauss_data")
  data <- WCCCGauss_data
}

## ----echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
 # Append dummies
 set.seed(789)
 update_data <- function(data_lst, num_dummies = p) {
   p <- ncol(data_lst$X)
   n <- nrow(data_lst$X)
 
   data_lst$X <- cbind(data_lst$X,
     matrix(complex(real = stats::rnorm(n * num_dummies, sd = 1),
                    imag = stats::rnorm(n * num_dummies, sd = 1)),
            nrow = n, ncol = num_dummies)
   )
   
   data_lst$num_dummies <- num_dummies
   data_lst$n <- n
   data_lst$p <- p
 
   return(data_lst)
 }
 data <- update_data(data, num_dummies = ncol(data$X))

## ----echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
# Run complex terminating lars algorithm
# ----------------------------------------
ctlars_obj <- ctlars::ctlars$new(data$X,
                                  data$y,
                                  num_dummies = data$num_dummies
                                  )

## ----terminated_solution_t_1, echo=TRUE, fig.align='center', message=TRUE, fig.width = 10, fig.height = 6, out.width = "90%"----
ctlars_obj$execute_clars_step(t_stop = 1, early_stop = TRUE, use_chol = TRUE)

## ----echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
# Test
# ----------------------------------------
cat("Proposed index set by CLARS: ", ctlars_obj$get_active_set(), "\n")
cat("True support: ", drop(data$support), "\n")
cat("Proposal in true support: ",
     intersect(ctlars_obj$get_active_set(), drop(data$support)), "\n")
cat("Dummies in proposal: ", ctlars_obj$get_active_dummies(), "\n")
 
plot(ctlars_obj)

## ----Warm_Start, echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
ctlars_obj$execute_clars_step(t_stop = 5)
plot(ctlars_obj)

## ----terminated_solution_t_5, echo=TRUE, fig.align='center', message=TRUE, fig.width=10, fig.height=6, out.width="90%"----
# Try a second example with t_stop = 5
# ---------------------------------------
ctlars_obj <- ctlars$new(data$X,
                         data$y,
                         has_intercept = FALSE,
                         num_dummies = data$num_dummies)

ctlars_obj$execute_clars_step(t_stop = 5, early_stop = TRUE)
plot(ctlars_obj)

## ----eval = TRUE--------------------------------------------------------------
ctlars_obj <- ctlars$new(data$X,
                          data$y,
                          has_intercept = FALSE,
                          num_dummies = data$num_dummies, 
                          verbose = FALSE)
 
ctlars_obj$execute_clars_step(t_stop = 1, early_stop = FALSE)

plot(ctlars_obj)

## -----------------------------------------------------------------------------
 # Evaluate FDP and TPP
 # ----------------------------------------
 chosen <- ctlars_obj$get_active_set()
 # remove dummies
 chosen <- setdiff(
   chosen, ctlars_obj$get_active_dummies()
 )
 FDP(chosen, drop(data$support))
 TPP(chosen, drop(data$support))
 
 # Plot R2 and SSR
 # ------------------
 plot(ctlars_obj$get_ssr(), type = "l", xlab = "Lars Step", ylab = "SSR")
 grid()
 
 plot(ctlars_obj$get_r2(), type = "l", xlab = "Lars Step", ylab = "R2")
 grid()

## ----FDR_and_TPR_sim, eval=TRUE, echo=TRUE, fig.align='center', message=FALSE, fig.width = 10, fig.height = 5, out.width = "90%"----
library(ctlars)

t_vec <- c(1, 2, 5, 10, 20, 50, 100)
sim_res <- simulate_model_TPR_FDR(
   t_vec,
   MC = 500,
   n = 100,
   p = 300,
   num_actives = 10,
   num_dummies = 300,
   set_snr = FALSE,
   verbose = FALSE
   )


## ----FDR_and_TPR_plot, eval=TRUE, echo=FALSE, fig.align='center', message=FALSE, fig.width = 10, fig.height = 5, out.width = "90%"----
# Plot results
library(ggplot2)
library(patchwork)

# data frame containing data to be plotted (FDR and TPR in %)
# -------------------------------------------------------------
plot_data <- data.frame(T_vec = t_vec,
                        FDR = 100 * sim_res$FDR,
                        TPR = 100 * sim_res$TPR
                        )
                        
 
# FDR vs. T
# ------------
fdr_vs_t <- ggplot(plot_data, aes(x = t_vec, y = FDR)) +
   labs(x = "T", y = "FDR") +
   scale_x_continuous(breaks = t_vec[-2], minor_breaks = c(2),
                      limits = c(t_vec[1], t_vec[length(t_vec)])) +
   scale_y_continuous(breaks = seq(0, 100, by = 10),
                      minor_breaks = c(), limits = c(0, 100)) +
   geom_line(linewidth = 1.5, colour = "#336C68") +
   geom_point(size = 2.5, colour = "#336C68") +
   theme_bw(base_size = 16) +
   theme(
     panel.background = element_rect(fill = "white",
                                     color = "black",
                                     linewidth = 1)
   ) +
   coord_fixed(ratio =  0.85 * t_vec[length(t_vec)] / (100 - 0))


# TPR vs. T
# ------------
tpr_vs_t <- ggplot(plot_data, aes(x = t_vec, y = TPR)) +
   labs(x = "T", y = "TPR") +
   scale_x_continuous(breaks = t_vec[-2], minor_breaks = c(2),
                      limits = c(t_vec[1], t_vec[length(t_vec)])) +
   scale_y_continuous(breaks = seq(0, 100, by = 10),
                      minor_breaks = c(), limits = c(0, 100)) +
   geom_line(linewidth = 1.5, colour = "#336C68") +
   geom_point(size = 2.5, colour = "#336C68") +
   theme_bw(base_size = 16) +
   theme(
     panel.background = element_rect(fill = "white",
                                     color = "black",
                                     linewidth = 1)
   ) +
   coord_fixed(ratio =  0.85 * t_vec[length(t_vec)] / (100 - 0))
 
fdr_vs_t + tpr_vs_t
 
 
# TPR vs. FDR
# ------------
tpr_vs_fdr <- ggplot(plot_data, aes(x = FDR, y = TPR)) +
   labs(x = "FDR",
        y = "TPR") +
   geom_line(linewidth = 1.5, colour = "#336C68") +
   geom_point(size = 2.5, colour = "#336C68") +
   theme_bw(base_size = 16) +
   theme(panel.background = element_rect(fill = "white",
                                         color = "black",
                                         linewidth = 1))
tpr_vs_fdr
 

