# Regression check for the two-start fit in PBT_var2_calc_MLE
# Uses made-up strata only. Each stratum has tagged (PBT-assigned) AI fish by release group and
# genetic stock, plus Unassigned AI fish by stock, with a tag rate of 0.9456 for every group.
# For each stratum the wild (Unassigned row) stock proportions from PBT_var2_calc_MLE are
# compared with the accounting solution (Unassigned counts minus expected untagged fish).
# Before the patch, 14 of these 32 strata returned a wild split pushed to the boundary
# (a stock at ~0% or ~100%) while reporting convergence.

library(escapeLGD)
calc <- escapeLGD:::PBT_var2_calc_MLE

stocks <- c("FALL", "HELLSC", "MFSALM", "SFSALM", "UPSALM")
full <- function(x) { v <- setNames(rep(0, length(stocks)), stocks); v[names(x)] <- x; v }

check_stratum <- function(hat, un, label, tag = 0.9456) {
	m <- rbind(do.call(rbind, lapply(hat, full)), Unassigned = full(un))
	rownames(m)[seq_along(hat)] <- paste0("G", seq_along(hat))
	n_hat <- rowSums(m)[-nrow(m)]
	untag <- n_hat * (1 - tag) / tag
	piGroup <- c(n_hat / tag, sum(un) - sum(untag)); piGroup <- piGroup / sum(piGroup)
	tagRates <- c(rep(tag, length(hat)), 1)
	wild <- calc(m, piGroup, tagRates)["Unassigned", ]
	sub <- colSums(sweep(m[-nrow(m), , drop = FALSE], 1, untag / n_hat, "*"))
	acct <- pmax(m["Unassigned", ] - sub, 0); acct <- acct / sum(acct)
	data.frame(stratum = label, tagged = sum(n_hat), unassigned = sum(un),
				  MLE_FALL = round(100 * wild["FALL"], 1), accounting_FALL = round(100 * acct["FALL"], 1),
				  max_diff_points = round(100 * max(abs(wild - acct)), 2), row.names = NULL)
}

un10 <- c(FALL = 55, HELLSC = 14, MFSALM = 7, SFSALM = 13, UPSALM = 6)
un7  <- c(FALL = 3, HELLSC = 40, MFSALM = 11, SFSALM = 31, UPSALM = 15)
unH  <- c(FALL = 6, HELLSC = 55, MFSALM = 7, SFSALM = 13, UPSALM = 14)
one <- function(n, k, s) lapply(rep(n %/% k, k), function(x) setNames(x, s))
mixed <- function(n) list(c(SFSALM = round(n*.4)), c(UPSALM = round(n*.2)), c(HELLSC = round(n*.2)),
								  c(MFSALM = round(n*.1)), c(FALL = round(n*.1)))
tiny <- list(c(SFSALM = 2), c(UPSALM = 1), c(HELLSC = 1), c(MFSALM = 1), c(UPSALM = 1), c(SFSALM = 1))

res <- rbind(
	do.call(rbind, lapply(c(1, 3, 10, 30), function(k) do.call(rbind, lapply(c(30, 60, 90, 120, 150, 300),
		function(n) check_stratum(one(n, k, "FALL"), un10, sprintf("%d FALL groups", k)))))),
	check_stratum(mixed(40), un10, "mixed-stock groups, FALL-heavy Unassigned"),
	check_stratum(mixed(146), un7, "mixed-stock groups, HELLSC/SFSALM Unassigned"),
	check_stratum(c(one(140, 10, "FALL"), tiny), un10, "140 FALL + 7 small non-FALL groups"),
	check_stratum(one(140, 10, "FALL"), un10, "140 FALL, no small groups"),
	do.call(rbind, lapply(c(40, 90, 140, 200), function(n) check_stratum(one(n, 10, "HELLSC"), unH, "10 HELLSC groups")))
)
print(res)
cat("\nStrata where MLE and accounting wild splits differ by more than 0.5 percentage points:",
	 sum(res$max_diff_points > 0.5), "of", nrow(res), "\n")
stopifnot(all(res$max_diff_points <= 0.5))
