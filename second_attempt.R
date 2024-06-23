rm(list = ls())

library(DAISIE)

# Load the data
data(frogs_datatable)

# Island parameters
island_age <- 30
mainland_pool <- 100L

# Prepare the data for DAISIE
datalist <- DAISIE_dataprep(frogs_datatable, island_age, M = mainland_pool)

# Parameter values
lambda_c <- 0.18
mu <- 0.03
K <- 0
gamma <- 0.0006
lambda_a <- 2

# How many species making it to the present?
k <- 0L

# How many unobserved species are allowed at most?
max_n_unobserved <- 100L

# Possible numbers of unobserved species
n <- 0:max_n_unobserved

# Prepare padding for boundary conditions
left_padding <- 2L
right_padding <- 1L

# Total number of state variables (incl. padding)
n_variables <- left_padding + length(n) + right_padding

# Set-up initial vectors of state variables
Qkn <- rep(0, n_variables)
QMkn <- rep(0, n_variables)


# Where does the system truly start?
id_n0 <- left_padding + 1L

# What is the range we are really concerned with?
n

# There is no observed nor unobserved species at the start
Qkn[id_n0] <- 1



dQkn <- mu * QMkn + lambda_a * QMkn[i - 1] + lambda_c * QMkn[i - 2] +
  lambda_c * (n + 2 * k - 1) * Qkn[i - 1] + mu * (n + 1) * Qkn[i + 1] -
  (mu + lambda_c) * (n + k) * Qkn - gamma * Qkn

dQMkn <- gamma * Qkn + lambda_c * (n + 2 * k - 1) * QMkn[i - 1] +
  mu * (n + 1) * QMkn[i + 1] - (mu + lambda_c) * (n + k) * QMkn -
  (mu + lambda_a + lambda_c) * QMkn

# Now vectorize
