# wrappers for HNC expand and variable GSI

#' treats GSI as known
#' @export
HNC_expand <- function(trap, stratAssign_comp, boots = 2000,
							  pbt_var, timestep_var, physTag_var,
							  adclip_var, tagRates,
							  H_vars, HNC_vars, W_vars, wc_binom, method = c("Account", "MLE"), testing = c("manual", "complete", "none")){
	method <- match.arg(method)
	testing <- match.arg(testing)
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
		if(method == "Account"){
			stratRes <- HNC_expand_one_strat(trap = trap %>% filter(stratum == s), H_vars = H_vars,
														HNC_vars = HNC_vars,
														W_vars = W_vars,
														wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
														pbt_var = pbt_var, tagRates = tagRates)
		} else if (method == "MLE"){
			stratRes <- HNC_expand_one_strat_MLE(trap = trap %>% filter(stratum == s), H_vars = H_vars,
														HNC_vars = HNC_vars,
														W_vars = W_vars,
														wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
														pbt_var = pbt_var, tagRates = tagRates)
		} else {
			stop("Error in method selection")
		}
		if(is.null(stratRes)) stop("Need complete cases for all quantities in stratum ", s)
		full_breakdown <- full_breakdown %>% bind_rows(stratRes$estimates %>% mutate(stratum = s))
	}

	if(boots < 1) return(list(full_breakdown, NULL))

	boot_breakdown <- tibble()
	# bootstraps
	skipped <- 0
	totalIter <- 0
	for(b in 1:boots){
		for(s in unique(stratAssign_comp$stratum)){
			while(TRUE){
				# resample genotyped and not genotyped separately b/c genotyping is only done on a specified subsample
				bootData <- bind_rows(trap %>% filter(stratum == s, is.na(pbtAssign)) %>%
											 	sample_n(nrow(.), replace = TRUE),
											 trap %>% filter(stratum == s, !is.na(pbtAssign)) %>%
											 	sample_n(nrow(.), replace = TRUE)
				)
				if(method == "Account"){
					bootRes <- HNC_expand_one_strat(trap = bootData, H_vars = H_vars,
															  HNC_vars = HNC_vars,
															  W_vars = W_vars,
															  wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
															  pbt_var = pbt_var, tagRates = tagRates)
				} else if (method == "MLE"){
					bootRes <- HNC_expand_one_strat_MLE(trap = bootData, H_vars = H_vars,
																	HNC_vars = HNC_vars,
																	W_vars = W_vars,
																	wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
																	pbt_var = pbt_var, tagRates = tagRates)
				} else {
					stop("Error in method selection")
				}
				totalIter <- totalIter + 1
				if(totalIter %% 100 == 0 && skipped / totalIter > .1) message((skipped / totalIter) * 100,
					" % of iterations are being skipped due to missing complete cases. ",
					"You may want to interrupt this function and prune observations with missing data")
				if(!is.null(bootRes)) break
				skipped <- skipped + 1
			}

			# make sure groups not sampled are estimated at zero
			if(testing == "manual") bootRes$estimates <- full_breakdown %>% filter(stratum == s) %>%
					select(rear, var1, var2) %>%
					left_join(bootRes$estimates, by = c("rear", "var1", "var2")) %>%
				mutate(total = replace_na(total, 0))

			boot_breakdown <- boot_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = b))
		}
	}
	if(testing == "complete") boot_breakdown <- boot_breakdown %>% group_by(rear, stratum, var1, var2) %>%
		complete(rear, stratum, var1, var2, boot, fill = list(total = 0))

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
							  GSI_draws, n_point = 100, GSI_var, method = c("Account", "MLE")){
	method <- match.arg(method)
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
			if(method == "Account"){
				bootRes_point <- HNC_expand_one_strat(trap = bootData_point, H_vars = H_vars,
																  HNC_vars = HNC_vars,
																  W_vars = W_vars,
																  wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
																  pbt_var = pbt_var, tagRates = tagRates)
			} else if (method == "MLE"){
				bootRes_point <- HNC_expand_one_strat_MLE(trap = bootData_point, H_vars = H_vars,
																  HNC_vars = HNC_vars,
																  W_vars = W_vars,
																  wc_expanded = wc_binom[[1]] %>% filter(stratum == s) %>% pull(wc),
																  pbt_var = pbt_var, tagRates = tagRates)
			} else {
				stop("Error in method selection")
			}
			if(is.null(bootRes_point)) stop("Need complete cases for all quantities in stratum ", s)
			full_breakdown <- full_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = i))
		}
	}
	full_breakdown <- full_breakdown %>% group_by(rear, stratum, var1, var2) %>% summarise(total = mean(total))

	if(boots < 1) return(list(full_breakdown, NULL))

	boot_breakdown <- tibble()
	# bootstraps
	for(b in 1:boots){
		for(s in unique(stratAssign_comp$stratum)){
			while(TRUE){
				# resample genotyped and not genotyped separately b/c genotyping is only done on a specified subsample
				bootData <- bind_rows(trap %>% filter(stratum == s, is.na(pbtAssign)) %>%
											 	sample_n(nrow(.), replace = TRUE),
											 trap %>% filter(stratum == s, !is.na(pbtAssign)) %>%
											 	sample_n(nrow(.), replace = TRUE)
				)
				# update GSI assignments
				ind_matches <- match(GSI_draws$Ind, bootData$Ind)
				bootData[ind_matches,GSI_var] <- GSI_draws[[b+1]] # first column is ind names
				if(method == "Account"){
					bootRes <- HNC_expand_one_strat(trap = bootData, H_vars = H_vars,
															  HNC_vars = HNC_vars,
															  W_vars = W_vars,
															  wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
															  pbt_var = pbt_var, tagRates = tagRates)
				} else if (method == "MLE"){
					bootRes <- HNC_expand_one_strat_MLE(trap = bootData, H_vars = H_vars,
																	HNC_vars = HNC_vars,
																	W_vars = W_vars,
																	wc_expanded = wc_binom[[2]][b, which(wc_binom[[1]]$stratum == s)],
																	pbt_var = pbt_var, tagRates = tagRates)
				} else {
					stop("Error in method selection")
				}
				totalIter <- totalIter + 1
				if(totalIter %% 100 == 0 && skipped / totalIter > .1) message((skipped / totalIter) * 100,
					" % of iterations are being skipped due to missing complete cases. ",
					"You may want to interrupt this function and prune observations with missing data")
				if(!is.null(bootRes)) break
				skipped <- skipped + 1
			}

			boot_breakdown <- boot_breakdown %>% bind_rows(bootRes$estimates %>% mutate(stratum = s, boot = b))
		}
	}

	return(list(full_breakdown, boot_breakdown))
}

