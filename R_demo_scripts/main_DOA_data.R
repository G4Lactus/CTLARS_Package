# Parameters
M <- 100    # Number of sensors
N <- 180    # Number of potential directions
K <- 5      # True number of sources
SNR <- 10   # Signal-to-Noise Ratio in dB

# Generate data
# True DOAs
theta_true <- sort(sample(N, K))
# Steering matrix
A <- exp(1i * 2 * pi * outer(0:(M-1), sin(pi * (0:(N-1)) / 180)))
# True source signals (single snapshot)
S_true <- (rnorm(K) + 1i * rnorm(K)) / sqrt(2)
# Received signal vector
X <- A[, theta_true] %*% S_true

# Add noise
noise_power <- 10^(-SNR/10)
X <- X + sqrt(noise_power/2) * (rnorm(M) + 1i * rnorm(M))  # Single snapshot

num_dummies <- ncol(A)
n <- nrow(A)
A_exp <- cbind(A,
               matrix(#complex(real = stats::rnorm(n * num_dummies, sd = 1),
                      #        imag = stats::rnorm(n * num_dummies, sd = 1)),
                      # circle border
                      exp(1i * stats::rnorm(n*num_dummies, 0, 2 * pi)),
                      nrow = n, ncol = num_dummies) )

# Run complex terminating lars algorithm
# ----------------------------------------
ctlars_obj <- ctlars$new(A_exp,
                         X,
                         has_intercept = TRUE,
                         standardize = TRUE,
                         num_dummies = num_dummies,
                         verbose = TRUE)

ctlars_obj$execute_clars_step(t_stop = 4)

cat("True actives: ", theta_true, "\n")
cat("Active set: ", ctlars_obj$get_active_set(), "\n")
cat("Intersection: ", intersect(ctlars_obj$get_active_set(), theta_true))
