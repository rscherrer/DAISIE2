# First set of derivatives
dQkn <- with(pars,
             mu * QMkn[ii] + lambda_a * QMkn[ii - 1] + lambda_c * QMkn[ii - 2] +
               lambda_c * (n + 2 * k - 1) * Qkn[ii - 1] + mu * (n + 1) * Qkn[ii + 1] -
               (mu + lambda_c) * (n + k) * Qkn[ii] - gamma * Qkn[ii]
)

# Second set of derivatives
dQMkn <- with(pars,
              gamma * Qkn[ii] + lambda_c * (n + 2 * k - 1) * QMkn[ii - 1] +
                mu * (n + 1) * QMkn[ii + 1] - (mu + lambda_c) * (n + k) * QMkn[ii] -
                (mu + lambda_a + lambda_c) * QMkn[ii]
)
