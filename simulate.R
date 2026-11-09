# Run with Rscript --vanilla simulate.R from the repository root.
# All observations are generated here. No input data or packages are needed.
RNGversion("4.0.0")
set.seed(20261109)

n <- 320L
c0 <- rnorm(n)
c0 <- c0 - mean(c0)                  # one vote per person
a <- 0.35 * c0 + rnorm(n)
m <- 0.30 + 0.50 * a + 0.40 * c0 + rnorm(n, sd = 0.70)

# A correlated random intercept and slope, independent of a, m and c0.
z0 <- rnorm(n)
z1 <- rnorm(n)
rho <- 0.30
intercept <- 2 + 0.70 * z0
slope <- 0.10 + 0.30 * a + 0.25 * m + 0.20 * a * m + 0.15 * c0 +
  0.25 * (rho * z0 + sqrt(1 - rho^2) * z1)

n_visits <- sample(3:6, n, replace = TRUE)
id <- rep(seq_len(n), n_visits)
occasion <- sequence(n_visits) - 1L
time <- occasion + ifelse(occasion == 0, 0, runif(length(id), -0.2, 0.2))
y <- intercept[id] + slope[id] * time + rnorm(length(id), sd = 0.45)

visits <- data.frame(id, time, y, a = a[id], m = m[id], c = c0[id])
# Fixed decimal precision and a headerless, space-delimited Mplus input file.
lines <- with(visits, sprintf("%d %.8f %.8f %.8f %.8f %.8f", id, time, y, a, m, c))
writeLines(lines, "synthetic.dat", useBytes = TRUE)
cat(sprintf("Generated synthetic.dat with %d people and %d visits.\n", n, nrow(visits)))
