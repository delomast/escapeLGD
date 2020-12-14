# now estimate composition
# this is just going to be the function for no expanding HNC into W and no using uncertainty in GSI
#' @export
ascension_composition <- function(trap, stratAssign_comp, boots = 2000,
											 pbt_var = NULL, tagRates = NULL,
											 H_vars, HNC_vars, W_vars, wc_binom){

	# need check for colomn names in trap (can't include "stratum"), or special handling if it does

	trap <- trap %>% left_join(stratAssign_comp, by = "sWeek")
	allStrat <- unique(c(stratAssign_comp$stratum))

	AD_rates <- list(tibble())
	AD_rates[[2]] <- matrix(nrow = boots, ncol = length(allStrat))
	AIphystag_rates <- list(tibble())
	AIphystag_rates[[2]] <- matrix(nrow = boots, ncol = length(allStrat))
	pbtOnly_rates <- list(tibble())
	pbtOnly_rates[[2]] <- matrix(nrow = boots, ncol = length(allStrat))
	seed_for_HNC_W <- c()
	for(i in 1:length(allStrat)){ # need to make sure all strata have values for all parameters
		s <- allStrat[i]
		tempTrap <- trap %>% filter(stratum == s, LGDMarkAD %in% c("AD", "AI"))
		if(nrow(tempTrap) == 0) stop("No valid trapped fish in stratum ", s)
		# AD vs AI
		AD_rates[[1]] <- AD_rates[[1]] %>% bind_rows(tibble(stratum = s,
																			 pClip = sum(tempTrap$LGDMarkAD == "AD") / nrow(tempTrap),
																			 totalTrap = nrow(tempTrap)))
		AD_rates[[2]][,i] <- rbinom(boots, AD_rates[[1]]$totalTrap[i], AD_rates[[1]]$pClip[i]) / AD_rates[[1]]$totalTrap[i]

		# split AI into phystag and no phystag
		tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI")
		if(nrow(tempTrap) == 0){
			AIphystag_rates[[1]] <- AIphystag_rates[[1]] %>%
				bind_rows(tibble(stratum = s, pPhys = 0, totalTrap = 0))
			AIphystag_rates[[2]][,i] <- 0
		} else {
			AIphystag_rates[[1]] <- AIphystag_rates[[1]] %>%
				bind_rows(tibble(stratum = s, pPhys = sum(tempTrap$physTag) / nrow(tempTrap),
									  totalTrap = nrow(tempTrap)))
			AIphystag_rates[[2]][,i] <- rbinom(boots, AIphystag_rates[[1]]$totalTrap[i],
														  AIphystag_rates[[1]]$pPhys[i]) / AIphystag_rates[[1]]$totalTrap[i]
		}

		# split no phystag into PBT-only HNC and unmarked, untagged
		tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI", !physTag)
		if(nrow(tempTrap) == 0){
			pbtOnly_rates[[1]] <- pbtOnly_rates[[1]] %>%
				bind_rows(tibble(stratum = s, pPBTonly = 0, totalTrap = 0))
			pbtOnly_rates[[2]][,i] <- 0
		} else {
			if(!any(!is.na(tempTrap$pbtAssign))){
				warning("all AI non-physTag samples in strata ", s, " were not assessed by PBT. Assuming all are Unassigned.")
				tempTrap$pbtAssign <- FALSE
			}
			tempTrap <- tempTrap %>% filter(!is.na(pbtAssign))
			# no PBT expansion
			pbtOnly_rates[[1]] <- pbtOnly_rates[[1]] %>%
				bind_rows(tibble(stratum = s, pPBTonly = sum(tempTrap$pbtAssign) / nrow(tempTrap),
									  totalTrap = nrow(tempTrap)))
			pbtOnly_rates[[2]][,i] <- rbinom(boots, pbtOnly_rates[[1]]$totalTrap[i],
														pbtOnly_rates[[1]]$pPBTonly[i]) / pbtOnly_rates[[1]]$totalTrap[i]

		}
	}

	# now multiply to get H, HNC, W breakdown
	rear_rates <- list() # point est, then H (cols are strata, rows boots), HNC, W
	rear_rates[[1]] <- tibble(stratum = AD_rates[[1]]$stratum,
									  H = AD_rates[[1]]$pClip,
									  HNC = (1 - AD_rates[[1]]$pClip) * # not clipped and (phystag or pbtOnly)
									  	(AIphystag_rates[[1]]$pPhys + ((1 - AIphystag_rates[[1]]$pPhys) * pbtOnly_rates[[1]]$pPBTonly))) %>%
		mutate(W = 1 - H - HNC)
	rear_rates[[2]] <- matrix(nrow = boots, ncol = nrow(rear_rates[[1]])) # bootstraps by strata for H
	rear_rates[[3]] <- matrix(nrow = boots, ncol = nrow(rear_rates[[1]])) # bootstraps by strata for HNC
	rear_rates[[4]] <- matrix(nrow = boots, ncol = nrow(rear_rates[[1]])) # bootstraps by strata for W
	for(i in 1:nrow(rear_rates[[1]])){
		rear_rates[[2]][,i] <- rbinom(boots, AD_rates[[1]]$totalTrap[i], AD_rates[[1]]$pClip[i]) / AD_rates[[1]]$totalTrap[i]
		temp_phystag <- rbinom(boots, AIphystag_rates[[1]]$totalTrap[i], AIphystag_rates[[1]]$pPhys[i]) /
			AIphystag_rates[[1]]$totalTrap[i]
		rear_rates[[3]][,i] <- (1 - rear_rates[[2]][,i]) * # not clipped and (phystag or pbtOnly)
			(temp_phystag + ((1 - temp_phystag) * (rbinom(boots, pbtOnly_rates[[1]]$totalTrap[i], pbtOnly_rates[[1]]$pPBTonly[i]) /
																	pbtOnly_rates[[1]]$totalTrap[i])))
		rear_rates[[4]][,i] <- 1 - rear_rates[[2]][,i] - rear_rates[[3]][,i]
	}

	# now breaking down rear types by given variables

	# first index is stratum, second is variable, then in third index first entry is
	# point estimates, rest are bootstraps
	# but for var2, third index is subgroup of var1, then fourth index is point and bootstraps
	H_break <- list()
	HNC_break <- list()
	W_break <- list()
	if(length(H_vars) < 1 || length(H_vars) > 2) stop("H_vars must be one or two column names")
	if(length(HNC_vars) < 1 || length(H_vars) > 2) stop("HNC_vars must be one or two column names")
	if(length(W_vars) < 1 || length(H_vars) > 2) stop("W_vars must be one or two column names")
	for(i in 1:nrow(rear_rates[[1]])){ # for each stratum
		s <- rear_rates[[1]]$stratum[i]
		tempTrap <- trap %>% filter(stratum == s, LGDMarkAD %in% c("AD", "AI"))
		# H subgroup proportions
		H_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[tempTrap$LGDMarkAD == "AD",],
													  vars = H_vars, pbt_var = pbt_var, tagRates = tagRates, boots = boots)
		tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI")

		# don't adjust for tag rates between groups
		# HNC subgroup proportions
		HNC_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[(tempTrap$physTag |
																							  	(!is.na(tempTrap$pbtAssign) & tempTrap$pbtAssign)),],
														 vars = HNC_vars, pbt_var = pbt_var, tagRates = tagRates,
														 boots = boots)
		tempTrap <- tempTrap %>% filter(!physTag)
		if(!any(!is.na(tempTrap$pbtAssign))) tempTrap$pbtAssign <- FALSE # don't need to warn this time, already did above
		# W subgroup proportions
		W_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[!is.na(tempTrap$pbtAssign) & !tempTrap$pbtAssign,],
													  vars = W_vars, pbt_var = pbt_var, tagRates = tagRates,
													  boots = boots)
	}

	# multiply through to get breakdown of ascensions

	full_breakdown <- tibble()
	boot_breakdown <- tibble()
	for(i in 1:nrow(rear_rates[[1]])){ # for each stratum
		# var1
		# H
		all_fish <- wc_binom[[1]]$wc[i] * rear_rates[[1]]$H[i]
		full_breakdown <- full_breakdown %>% bind_rows(
			H_break[[i]][[1]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
														 rear = "H", var1 = group, var2 = NA) %>%
				select(stratum, rear, var1, var2, total)
		)
		boot_totals <- H_break[[i]][[1]][[2]] * (wc_binom[[2]][,i] * rear_rates[[2]][,i])
		# paranoid about element wise multiplication
		if(!all.equal(rowSums(boot_totals), (wc_binom[[2]][,i] * rear_rates[[2]][,i]))) warning("possible internal error")
		for(j in 1:ncol(boot_totals)){
			boot_breakdown <- bind_rows(boot_breakdown,
												 tibble(total = boot_totals[,j], var1 = colnames(boot_totals)[j],
												 		 stratum = rear_rates[[1]]$stratum[i],
												 		 var2 = NA, rear = "H", boot = 1:nrow(boot_totals)) %>%
												 	select(boot, stratum, rear, var1, var2, total))
		}

		# HNC
		all_fish <- wc_binom[[1]]$wc[i] * rear_rates[[1]]$HNC[i]
		full_breakdown <- full_breakdown %>% bind_rows(
			HNC_break[[i]][[1]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
															rear = "HNC", var1 = group, var2 = NA) %>%
				select(stratum, rear, var1, var2, total)
		)
		boot_totals <- HNC_break[[i]][[1]][[2]] * (wc_binom[[2]][,i] * rear_rates[[3]][,i])
		# paranoid about element wise multiplication
		if(!all.equal(rowSums(boot_totals), (wc_binom[[2]][,i] * rear_rates[[3]][,i]))) warning("possible internal error")
		for(j in 1:ncol(boot_totals)){
			boot_breakdown <- bind_rows(boot_breakdown,
												 tibble(total = boot_totals[,j], var1 = colnames(boot_totals)[j],
												 		 stratum = rear_rates[[1]]$stratum[i],
												 		 var2 = NA, rear = "HNC", boot = 1:nrow(boot_totals)) %>%
												 	select(boot, stratum, rear, var1, var2, total))
		}

		# W
		all_fish <- wc_binom[[1]]$wc[i] * rear_rates[[1]]$W[i]
		full_breakdown <- full_breakdown %>% bind_rows(
			W_break[[i]][[1]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
														 rear = "W", var1 = group, var2 = NA) %>%
				select(stratum, rear, var1, var2, total)
		)
		boot_totals <- W_break[[i]][[1]][[2]] * (wc_binom[[2]][,i] * rear_rates[[4]][,i])
		# paranoid about element wise multiplication
		if(!all.equal(rowSums(boot_totals), (wc_binom[[2]][,i] * rear_rates[[4]][,i]))) warning("possible internal error")
		for(j in 1:ncol(boot_totals)){
			boot_breakdown <- bind_rows(boot_breakdown,
												 tibble(total = boot_totals[,j], var1 = colnames(boot_totals)[j],
												 		 stratum = rear_rates[[1]]$stratum[i],
												 		 var2 = NA, rear = "W", boot = 1:nrow(boot_totals)) %>%
												 	select(boot, stratum, rear, var1, var2, total))
		}

		# var2
		if(length(H_vars) > 1){
			cats <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "H") %>%
				pull(var1) %>% unique
			for(c in cats){
				all_fish <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "H",
																  var1 == c) %>% pull(total)
				full_breakdown <- full_breakdown %>% bind_rows(
					H_break[[i]][[2]][[c]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
																		rear = "H", var1 = c, var2 = group) %>%
						select(stratum, rear, var1, var2, total)
				)
				all_fish <- boot_breakdown %>% filter(rear == "H", stratum == rear_rates[[1]]$stratum[i], var1 == c, is.na(var2)) %>%
					arrange(boot) %>% pull(total)
				boot_totals <- H_break[[i]][[2]][[c]][[2]] * all_fish
				# paranoid about element wise multiplication
				if(!all.equal(rowSums(boot_totals), all_fish)) warning("possible internal error")
				for(j in 1:ncol(boot_totals)){
					boot_breakdown <- bind_rows(boot_breakdown,
														 tibble(total = boot_totals[,j], var1 = c,
														 		 stratum = rear_rates[[1]]$stratum[i],
														 		 var2 = colnames(boot_totals)[j], rear = "H", boot = 1:nrow(boot_totals)) %>%
														 	select(boot, stratum, rear, var1, var2, total))
				}
			}
		}
		if(length(HNC_vars) > 1){
			cats <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "HNC") %>%
				pull(var1) %>% unique
			for(c in cats){
				all_fish <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "HNC",
																  var1 == c) %>% pull(total)
				full_breakdown <- full_breakdown %>% bind_rows(
					HNC_break[[i]][[2]][[c]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
																		  rear = "HNC", var1 = c, var2 = group) %>%
						select(stratum, rear, var1, var2, total)
				)
				all_fish <- boot_breakdown %>% filter(rear == "HNC", stratum == rear_rates[[1]]$stratum[i], var1 == c, is.na(var2)) %>%
					arrange(boot) %>% pull(total)
				boot_totals <- HNC_break[[i]][[2]][[c]][[2]] * all_fish
				# paranoid about element wise multiplication
				if(!all.equal(rowSums(boot_totals), all_fish)) warning("possible internal error")
				for(j in 1:ncol(boot_totals)){
					boot_breakdown <- bind_rows(boot_breakdown,
														 tibble(total = boot_totals[,j], var1 = c,
														 		 stratum = rear_rates[[1]]$stratum[i],
														 		 var2 = colnames(boot_totals)[j], rear = "HNC", boot = 1:nrow(boot_totals)) %>%
														 	select(boot, stratum, rear, var1, var2, total))
				}
			}
		}
		if(length(W_vars) > 1){
			cats <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "W") %>%
				pull(var1) %>% unique
			for(c in cats){
				all_fish <- full_breakdown %>% filter(stratum == rear_rates[[1]]$stratum[i], rear == "W",
																  var1 == c) %>% pull(total)
				full_breakdown <- full_breakdown %>% bind_rows(
					W_break[[i]][[2]][[c]][[1]] %>% mutate(total = prop * all_fish, stratum = rear_rates[[1]]$stratum[i],
																		rear = "W", var1 = c, var2 = group) %>%
						select(stratum, rear, var1, var2, total)
				)
				all_fish <- boot_breakdown %>% filter(rear == "W", stratum == rear_rates[[1]]$stratum[i], var1 == c, is.na(var2)) %>%
					arrange(boot) %>% pull(total)
				boot_totals <- W_break[[i]][[2]][[c]][[2]] * all_fish
				# paranoid about element wise multiplication
				if(!all.equal(rowSums(boot_totals), all_fish)) warning("possible internal error")
				for(j in 1:ncol(boot_totals)){
					boot_breakdown <- bind_rows(boot_breakdown,
														 tibble(total = boot_totals[,j], var1 = c,
														 		 stratum = rear_rates[[1]]$stratum[i],
														 		 var2 = colnames(boot_totals)[j], rear = "W", boot = 1:nrow(boot_totals)) %>%
														 	select(boot, stratum, rear, var1, var2, total))
				}
			}
		}
	}

	return(list(full_breakdown, boot_breakdown))
}
