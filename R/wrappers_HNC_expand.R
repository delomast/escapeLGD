# wrappers for HNC expand and variable GSI

#' @export
HNC_expand <- function(trap, stratAssign_comp, boots = 2000,
							  pbt_var, tagRates,
							  H_vars, HNC_vars, W_vars, wc_binom){
	# make sure tagRates conforms
	colnames(tagRates) <- c("group", "tagRate")
	if("Unassigned" %in% tagRates[[1]]){
		if(tagRates[[2]][tagRates[[1]] == "Unassigned"] != 1){
			stop("Unassigned must not be included in the tag rate file or have a tag rate of 1")
		}
	} else {
		# add Unassigned with tag rate of 1
		tagRates <- bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1))
	}


	# checking inputs




	full_breakdown <- tibble()
	# point estimate
	for(s in unique(stratAssign_comp$stratum)){
		stratRes <- HNC_expand_one_strat(trap = trap %>% filter(stratum == s), H_vars = H_vars,
									HNC_vars = HNC_vars,
									W_vars = W_vars,
									wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
									pbt_var = pbt_var, tagRates = tagRates)
		full_breakdown <- full_breakdown %>% bind_rows(stratRes$estimates %>% mutate(stratum = s))
	}

	boot_breakdown <- tibble()
	# bootstraps
	for(b in 1:boots){
		for(s in unique(stratAssign_comp$stratum)){
			# resample genotyped and not genotyped separately b/c genotyping is only done on a specified subsample
			bootData <- bind_rows(trap %>% filter(stratum == s, is.na(pbtAssign)) %>%
				sample_n(nrow(.), replace = TRUE),
				trap %>% filter(stratum == s, !is.na(pbtAssign)) %>%
					sample_n(nrow(.), replace = TRUE)
				)
			bootRes <- HNC_expand_one_strat(trap = bootData, H_vars = H_vars,
														HNC_vars = HNC_vars,
														W_vars = W_vars,
														wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
														pbt_var = pbt_var, tagRates = tagRates)
			boot_breakdown <- boot_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = b))
		}
	}

	return(list(full_breakdown, boot_breakdown))
}
