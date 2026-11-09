# Optional base R summary for the supplied Mplus input, run after estimation.
output <- readLines("four_way.out", warn = FALSE)
if (!any(trimws(output) == "THE MODEL ESTIMATION TERMINATED NORMALLY")) {
  stop("Check four_way.out. Mplus did not report normal termination.")
}

# Match names in the output to columns in BPARAMETERS. The first two
# columns identify the chain and iteration. Model edits can shift the rest.
components <- c("CDE", "INTREF", "INTMED", "PIE")
effects <- c(components, "TE")
parameter_lines <- grep("^ *Parameter [0-9]+, ", output, value = TRUE)
parameter_names <- sub("^ *Parameter [0-9]+, ", "", parameter_lines)
parameter_numbers <- as.integer(sub("^ *Parameter ([0-9]+), .*", "\\1", parameter_lines))
columns <- parameter_numbers[match(effects, parameter_names)] + 2L
if (anyNA(columns)) stop("Missing a decomposition parameter in four_way.out.")

draws <- read.table("posterior.dat", header = FALSE)
if (ncol(draws) != max(parameter_numbers) + 2L) {
  stop("The posterior file does not match the parameter list in four_way.out.")
}
# For this fixed-iteration run Mplus retains the second half of each chain.
keep <- draws[[2]] > ave(draws[[2]], draws[[1]], FUN = max) / 2
values <- as.matrix(draws[keep, columns])
colnames(values) <- effects
if (!all(is.finite(values))) stop("Non-finite posterior values.")
parts <- values[, components, drop = FALSE]
if (any(abs(rowSums(parts) - values[, "TE"]) > 5e-5)) {
  stop("Components did not sum to TE within the precision of the saved draws.")
}

intervals <- function(x) {
  t(apply(x, 2, quantile, probs = c(0.5, 0.025, 0.975), names = FALSE))
}
effect_table <- intervals(values)
colnames(effect_table) <- c("Median", "Lower 95%", "Upper 95%")
cat("Difference in outcome units per time unit\n")
print(round(effect_table, 3))

# Magnitude shares sum to 100% per draw and ignore the component signs.
# An all-zero draw has undefined shares.
denominator <- rowSums(abs(parts))
defined <- denominator > 0
if (any(!defined)) message(sum(!defined), " draws have undefined magnitude shares.")
if (any(defined)) {
  shares <- 100 * abs(parts[defined, , drop = FALSE]) / denominator[defined]
  share_table <- intervals(shares)
  colnames(share_table) <- colnames(effect_table)
  cat("\nShare of the sum of absolute component magnitudes (%)\n")
  print(round(share_table, 1))
}
