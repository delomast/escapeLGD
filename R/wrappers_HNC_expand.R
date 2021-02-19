# wrappers for HNC expand and variable GSI

#' treats GSI as known
#' @export
HNC_expand <- function(trap, stratAssign_comp, boots = 2000,
							  pbt_var, timestep_var, physTag_var,
							  adclip_var, tagRates,
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
	if(length(H_vars) < 1 || length(H_vars) > 2) stop("H_vars must be one or two column names")
	if(length(HNC_vars) < 1 || length(HNC_vars) > 2) stop("HNC_vars must be one or two column names")
	if(length(W_vars) < 1 || length(W_vars) > 2) stop("W_vars must be one or two column names")
	allVars <- unique(c(H_vars, HNC_vars, W_vars))
	if(any(!allVars %in% colnames(trap))) stop("*_vars names must be column names in trap")
	if(any(c(timestep_var, pbt_var, physTag_var, allVars) == "LGDMarkAD")) stop("LGDMarkAD is reserved for adclip_var or must not be used")
	if(any(c(timestep_var, pbt_var, adclip_var, allVars) == "physTag")) stop("physTag is reserved for physTag_var or must not be used")
	if(any(c(physTag_var, pbt_var, adclip_var, allVars) == "sWeek")) stop("sWeek is reserved for timestep_var or must not be used")

	# check release groups in tag rate file

	# select relevant columns
	trap <- trap %>% select(all_of(c(timestep_var, pbt_var, adclip_var, physTag_var, allVars))) %>%
		rename(LGDMarkAD = !!as.symbol(adclip_var), physTag = !!as.symbol(physTag_var), sWeek = !!as.symbol(timestep_var))
	# make sure no reserved column names are used
	if(any(colnames(trap) %in% c("stratum", "pbtAssign"))) stop("stratum and pbtAssign are reserved column names and must not be used as column names in trap")

	# add "stratum" column to trap using stratAssign_comp and add pbtAssign
	colnames(stratAssign_comp)[1:2]  <- c("sWeek", "stratum")
	trap <- trap %>% mutate(pbtAssign = TRUE) %>% left_join(stratAssign_comp, by = "sWeek")
	trap$pbtAssign[trap[[pbt_var]] == "Unassigned"] <- FALSE
	trap$pbtAssign[is.na(trap[[pbt_var]])] <- NA


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
														wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
														pbt_var = pbt_var, tagRates = tagRates)
			boot_breakdown <- boot_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = b))
		}
	}

	return(list(full_breakdown, boot_breakdown))
}



#' Treats GSI as unknown
#' @param GSI_draws a tibble wiht the first column individuals names, the rest GSI
#'   values fro draws from the posterior (each column a draw). Used in order.
#' @param n_point number of iterations in the bootstrap loop (only resampling GSI values)
#'   for calculating point estimates
#' @param GSI_var the column name in \code{trap} that has the GSI values.
#' @export
HNC_expand_unkGSI <- function(trap, stratAssign_comp, boots = 2000,
							  pbt_var, tagRates,
							  H_vars, HNC_vars, W_vars, wc_binom,
							  GSI_draws, n_point = 100, GSI_var){
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
	if(length(H_vars) < 1 || length(H_vars) > 2) stop("H_vars must be one or two column names")
	if(length(HNC_vars) < 1 || length(HNC_vars) > 2) stop("HNC_vars must be one or two column names")
	if(length(W_vars) < 1 || length(W_vars) > 2) stop("W_vars must be one or two column names")
	allVars <- unique(c(H_vars, HNC_vars, W_vars))
	if(any(!allVars %in% colnames(trap))) stop("*_vars names must be column names in trap")
	if(any(c(timestep_var, pbt_var, physTag_var, allVars) == "LGDMarkAD")) stop("LGDMarkAD is reserved for adclip_var or must not be used")
	if(any(c(timestep_var, pbt_var, adclip_var, allVars) == "physTag")) stop("physTag is reserved for physTag_var or must not be used")
	if(any(c(physTag_var, pbt_var, adclip_var, allVars) == "sWeek")) stop("sWeek is reserved for timestep_var or must not be used")
	if(any(!trap[[sampID]] %in% GSI_draws[[1]])) stop("all sampID in trap must be in GSI_draws")



	full_breakdown <- tibble()
	# point estimate
	bootData_point <- trap
	ind_matches <- match(GSI_draws$Ind, bootData_point$Ind)
	for(i in 1:n_point){
		# update GSI assignments
		bootData_point[ind_matches,GSI_var] <- GSI_draws[[i+1]] # first column is ind names
		for(s in unique(stratAssign_comp$stratum)){
			bootRes_point <- HNC_expand_one_strat(trap = bootData_point, H_vars = H_vars,
													  HNC_vars = HNC_vars,
													  W_vars = W_vars,
													  wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
													  pbt_var = pbt_var, tagRates = tagRates)
			full_breakdown <- full_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = i))
		}
	}
	full_breakdown <- full_breakdown %>% group_by(rear, stratum, var1, var2) %>% summarise(total = mean(total))


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
			# update GSI assignments
			ind_matches <- match(GSI_draws$Ind, bootData$Ind)
			bootData[ind_matches,GSI_var] <- GSI_draws[[b+1]] # first column is ind names
			bootRes <- HNC_expand_one_strat(trap = bootData, H_vars = H_vars,
													  HNC_vars = HNC_vars,
													  W_vars = W_vars,
													  wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
													  pbt_var = pbt_var, tagRates = tagRates)
			boot_breakdown <- boot_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = b))
		}
	}

	return(list(full_breakdown, boot_breakdown))
}

