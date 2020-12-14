# night time passage, fallback and reascension functions

#'
#' estimate fallback, reascension, and nighttime passage
#' @param stratAssign_fallback tibble with sWeek, stockGroup, stratum showing what stratum each sWeek
#'   corresponds to for each stockGroup
#' @param stratAssign_night tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for nighttime passage
#'
#' @export
nightFall <- function(wc, full_reascend, full_night, stratAssign_fallback, stratAssign_night,
							 boots = 2000, full_spillway = NULL){
	# some input checking and maybe column name changes

	# note that this allows different stratification for fallback by stockGroup, nighttime passage, and SCOBI
	# only constraint is that strata changes in fallback by stockGroup must coincide with strata changes in SCOBI

	full_reascend <- full_reascend %>% left_join(stratAssign_fallback, by = c("sWeek", "stockGroup")) %>%
		group_by(stockGroup, stratum) %>% summarise(numReascend = sum(numReascend), totalPass = sum(totalPass))
	full_night <- full_night %>% left_join(stratAssign_night, by = "sWeek") %>% group_by(stratum) %>%
		summarise(nightPass = sum(nightPass), totalPass = sum(totalPass))
	if(!is.null(full_spillway)){
		full_spillway <- full_spillway %>% left_join(stratAssign_fallback, by = c("sWeek", "stockGroup")) %>%
			group_by(stockGroup, stratum) %>% summarise(totalFall = sum(totalFall), laterAscend = sum(laterAscend))
		fallback_data <- full_reascend %>% full_join(full_spillway, by = c("stockGroup", "stratum"))
	} else {
		fallback_data <- full_reascend %>% mutate(totalFall = 0, laterAscend = 0)
	}

	# make sure all strata are present and nonzero totals
	allStrat <- unique(c(stratAssign_fallback$stratum, stratAssign_night$stratum))
	for(g in unique(fallback_data$stockGroup)){
		if(any(!allStrat %in% fallback_data$stratum[fallback_data$stockGroup == g])) stop("not all strata are in the fallback and reascend dataset for stockGroup ", g)
	}
	if(any(!allStrat %in% full_night$stratum)) stop("not all strata are in the nighttime passage dataset ")
	if(any(fallback_data$totalPass == 0)) stop("total of zero PIT tags passing ladder in fallback data for at least one stratum")
	if(any(full_night$totalPass == 0)) stop("total of zero PIT tags passing ladder in nighttime passage data for at least one stratum")

	# first entry is point estimates
	# second entry is matrix with rows as bootstrap iterations and columns corresponding
	# to the rows of the first entry (strata and stock group)
	fallback_rates <- list(fallback_data %>% mutate(p_fa = NA))
	fallback_rates[[2]] <- matrix(nrow = boots, ncol = nrow(fallback_rates[[1]]))
	for(i in 1:nrow(fallback_rates[[1]])){
		# if no data for P(reascend|fallback) [eg no spillway data], then
		# assume P(reascend|fallback) = 1 [eg no fallback w/o reascension]
		if(fallback_rates[[1]]$totalFall[i] == 0){
			warning("Assuming no fallback withOUT reascension for ", "stratum ",
					  fallback_rates[[1]]$stratum[i], " stockGroup ", fallback_rates[[1]]$stockGroup[i])
			fallback_rates[[1]]$p_fa[i] <- fallback_rates[[1]]$numReascend[i] / fallback_rates[[1]]$totalPass[i]
			fallback_rates[[2]][,i] <- rbinom(boots, fallback_rates[[1]]$totalPass[i],
														 fallback_rates[[1]]$p_fa[i]) / fallback_rates[[1]]$totalPass[i]
		} else { # otherwise infer both
			opts <- optim(par = c(.1,.9), fn = optimllh, gr = gradient_fallback_log_likelihood,
							  dfr = fallback_rates[[1]]$laterAscend[i], df = fallback_rates[[1]]$totalFall[i],
							  dr = fallback_rates[[1]]$numReascend[i], dt = fallback_rates[[1]]$totalPass[i],
							  control = list(fnscale = -1), method = "L-BFGS-B", upper = 1 - 1e-7,
							  lower = 1e-7)
			if(opts$convergence != 0) warning("Convergence error inferring P(fallback)")
			fallback_rates[[1]]$p_fa[i] <- opts$par[1]

			# bootstrap by resampling observations within each dataset (ladder_filtered and spillway_count) separately
			# note that this gives crazy answers for some small sample sizes, like all bootstrap routines
			###
			# here we are nonparametrically resampling observations within each dataset,
			# but because of the way the data is this is the same as sampling a
			# binomial random variable
			boot_reas_data <- rbinom(boots, fallback_rates[[1]]$totalPass[i],
											 fallback_rates[[1]]$numReascend[i] / fallback_rates[[1]]$totalPass[i])
			boot_fallback_data <- rbinom(boots, fallback_rates[[1]]$totalFall[i],
												  fallback_rates[[1]]$laterAscend[i] / fallback_rates[[1]]$totalFall[i])
			convergeFail <- 0
			for(b in 1:boots){
				b_opts <- optim(par = c(.1,.9), fn = optimllh, gr = gradient_fallback_log_likelihood,
									 dfr = boot_fallback_data[b], df = fallback_rates[[1]]$totalFall[i],
									 dr = boot_reas_data[b], dt =  fallback_rates[[1]]$totalPass[i],
									 control = list(fnscale = -1), method = "L-BFGS-B", upper = 1 - 1e-7,
									 lower = 1e-7)
				if(b_opts$convergence != 0) convergeFail <- convergeFail + 1
				fallback_rates[[2]][b,i] <- b_opts$par[1]
			}
			if(convergeFail > 0) warning("A total of ", convergeFail, "convergence errors while bootstrapping ",
												  "stratum ", fallback_rates[[1]]$stratum[i], " stockGroup ",
												  fallback_rates[[1]]$stockGroup[i])
		}
	}


	# similar setup for nightime passage
	nightPassage_rates <- list(full_night %>% mutate(p_night = nightPass / totalPass))
	nightPassage_rates[[2]] <- matrix(nrow = boots, ncol = nrow(nightPassage_rates[[1]]))
	for(i in 1:nrow(nightPassage_rates[[1]])){
		nightPassage_rates[[2]][,i] <- rbinom(boots, nightPassage_rates[[1]]$totalPass[i],
														  nightPassage_rates[[1]]$p_night[i]) / nightPassage_rates[[1]]$totalPass[i]
	}

	return(list(fallback_rates = fallback_rates,
					nightPassage_rates = nightPassage_rates))
}


#' estimate and bootstrap window count as a binomial and expanding for nighttime passage
#' @param nightPassage_rates one of the outputs of \code{nightFall}
#' @param wc tibble with two columns sWeek and count of fish (NOT expanded for wc_prop)
#' @param wc_prop the proportion of the time fish are counted
#' @param stratAssign_night tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for nighttime passage
#' @param stratAssign_comp tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for composition estimation using the trap data
#' @export
expand_wc_binom_night <- function(nightPassage_rates, wc, wc_prop, stratAssign_comp, stratAssign_night,
											 boots = 2000){


	# first, expand each sWeek for wc_prop and nighttime passage

	# similar setup for window count
	wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop))) # expanding for wc_prop
	wc_binom[[2]] <- matrix(nrow = boots, ncol = nrow(wc_binom[[1]]))
	for(i in 1:nrow(wc_binom[[1]])){
		wc_binom[[2]][,i] <- rbinom(boots, wc_binom[[1]]$wc[i], wc_prop) / wc_prop
		# nighttime passage
		temp_nightStrata <- stratAssign_night$stratum[stratAssign_night$sWeek == wc$sWeek[i]]
		temp_nightStrata <- which(nightPassage_rates[[1]]$stratum == temp_nightStrata) # change into row/column number
		temp_nightRate <- nightPassage_rates[[1]]$p_night[temp_nightStrata]
		# expanding point estimate
		wc_binom[[1]]$wc[i] <-  wc_binom[[1]]$wc[i] / (1 - temp_nightRate)
		# expanding bootstrap
		wc_binom[[2]][,i] <- wc_binom[[2]][,i] / (1 - nightPassage_rates[[2]][,temp_nightStrata])
	}

	# then sum together according to strata
	wc_binom[[1]] <- wc_binom[[1]] %>% left_join(stratAssign_comp, by = "sWeek")

	temp_bootstrapMatrix <- matrix(0, nrow = boots, ncol = n_distinct(wc_binom[[1]]$stratum))
	pos <- 1
	for(i in 1:nrow(wc_binom[[1]])){
		if(i == 1){
			temp_bootstrapMatrix[,pos] <- temp_bootstrapMatrix[,pos] + wc_binom[[2]][,i]
		} else {
			if(wc_binom[[1]]$stratum[i-1] != wc_binom[[1]]$stratum[i]) pos <- pos + 1
			temp_bootstrapMatrix[,pos] <- temp_bootstrapMatrix[,pos] + wc_binom[[2]][,i]
		}
	}
	wc_binom[[1]] <- wc_binom[[1]] %>% group_by(stratum) %>% summarise(wc = sum(wc))
	wc_binom[[2]] <- temp_bootstrapMatrix

	return(wc_binom)
}

# now estimate composition
# this is just going to be the function for no expanding HNC into W and no using uncertainty in GSI
#' @export
ascension_composition <- function(trap, stratAssign_comp, boots = 2000,
											 pbt_var = NULL, tagRates = NULL,
											 H_vars, HNC_vars, W_vars){

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
		s <- rear_rates$stratum[i]
		tempTrap <- trap %>% filter(stratum == s, LGDMarkAD %in% c("AD", "AI"))
		# H subgroup proportions
		H_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[tempTrap$LGDMarkAD == "AD",],
													  vars = H_vars, pbt_var = pbt_var, tagRates = tagRates)
		tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI")

		# don't adjust for tag rates between groups
		# HNC subgroup proportions
		HNC_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[(tempTrap$physTag |
																							  	(!is.na(tempTrap$pbtAssign) & tempTrap$pbtAssign)),],
														 vars = HNC_vars, pbt_var = pbt_var, tagRates = tagRates)
		tempTrap <- tempTrap %>% filter(!physTag)
		if(!any(!is.na(tempTrap$pbtAssign))) tempTrap$pbtAssign <- FALSE # don't need to warn this time, already did above
		# W subgroup proportions
		W_break[[i]] <- subGroup_breakdown(trapStratumData = tempTrap[!is.na(tempTrap$pbtAssign) & !tempTrap$pbtAssign,],
													  vars = W_vars, pbt_var = pbt_var, tagRates = tagRates)
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

#' @export
apply_fallback_rates <- function(breakdown, fallback_rates,
											split_H_fallback = c("var1", "var2", "both"),
											split_HNC_fallback = c("var1", "var2", "both"),
											split_W_fallback = c("var1", "var2", "both"),
											H_groups = NULL, HNC_groups = NULL, W_groups = NULL,
											stratAssign_fallback, stratAssign_comp, alpha_ci = .1,
											output_type = c("summary", "W_boot", "full")){
	split_H_fallback <- match.arg(split_H_fallback)
	split_HNC_fallback <- match.arg(split_HNC_fallback)
	split_W_fallback <- match.arg(split_W_fallback)
	output_type <- match.arg(output_type)

	full_breakdown <- breakdown[[1]]
	boot_breakdown <- breakdown[[2]]
	rm(breakdown)

	group_tables <- list()
	if(split_H_fallback == "var1"){
		sc <- "var1"
	} else if (split_H_fallback == "var2"){
		sc <- "var2"
	} else {
		sc <- c("var1", "var2")
	}
	group_tables[[1]] <- full_breakdown %>% filter(rear == "H") %>%
		select(all_of(sc)) %>% distinct %>% mutate(stockGroup = NA)
	if("var2" %in% sc) group_tables[[1]] <- group_tables[[1]] %>% filter(!is.na(var2))
	if(split_HNC_fallback == "var1"){
		sc <- "var1"
	} else if (split_HNC_fallback == "var2"){
		sc <- "var2"
	} else {
		sc <- c("var1", "var2")
	}
	group_tables[[2]] <- full_breakdown %>% filter(rear == "HNC")  %>%
		select(all_of(sc)) %>% distinct %>% mutate(stockGroup = NA)
	if("var2" %in% sc) group_tables[[2]] <- group_tables[[2]] %>% filter(!is.na(var2))
	if(split_W_fallback == "var1"){
		sc <- "var1"
	} else if (split_W_fallback == "var2"){
		sc <- "var2"
	} else {
		sc <- c("var1", "var2")
	}
	group_tables[[3]] <- full_breakdown %>% filter(rear == "W")  %>%
		select(all_of(sc)) %>% distinct %>% mutate(stockGroup = NA)
	if("var2" %in% sc) group_tables[[3]] <- group_tables[[3]] %>% filter(!is.na(var2))
	names(group_tables) <- c("H", "HNC", "W")

	if(is.null(H_groups) || is.null(HNC_groups) || is.null(W_groups)){
		message("Returning a list of templates for assigning fallback rate groups")
		return(group_tables)
	}
	group_tables[[1]] <- group_tables[[1]] %>%
		left_join(H_groups, by = colnames(group_tables[[1]] %>% select(-stockGroup)))
	group_tables[[2]] <- group_tables[[2]] %>%
		left_join(HNC_groups, by = colnames(group_tables[[2]] %>% select(-stockGroup)))
	group_tables[[3]] <- group_tables[[3]] %>%
		left_join(W_groups, by = colnames(group_tables[[3]] %>% select(-stockGroup)))
	for(i in 1:3){
		if(any(is.na(group_tables[[i]]$stockGroup))){
			message("Missing an entry for ", c("H", "HNC", "W")[i], ". Returning a list of partially filled in templates.")
			return(group_tables)
		}
	}

	# matching up strata
	strataMatchUp <- stratAssign_comp %>% full_join(stratAssign_fallback %>% mutate(fallback = stratum), by = "sWeek")

	# now multiply fallback rates as appropriate

	full_breakdown_H <- full_breakdown %>% filter(rear == "H")
	full_breakdown_HNC <- full_breakdown %>% filter(rear == "HNC")
	full_breakdown_W <- full_breakdown %>% filter(rear == "W")
	rm(full_breakdown) # save memory
	boot_breakdown_H <- boot_breakdown %>% filter(rear == "H")
	boot_breakdown_HNC <- boot_breakdown %>% filter(rear == "HNC")
	boot_breakdown_W <- boot_breakdown %>% filter(rear == "W")
	rm(boot_breakdown) # save memory
	# removing unadjusted estimates
	if(split_H_fallback != "var1") {
		full_breakdown_H <- full_breakdown_H %>% filter(!is.na(var2))
		boot_breakdown_H <- boot_breakdown_H %>% filter(!is.na(var2))
	}
	if(split_HNC_fallback != "var1") {
		full_breakdown_HNC <- full_breakdown_HNC %>% filter(!is.na(var2))
		boot_breakdown_HNC <- boot_breakdown_HNC %>% filter(!is.na(var2))
	}
	if(split_W_fallback != "var1") {
		full_breakdown_W <- full_breakdown_W %>% filter(!is.na(var2))
		boot_breakdown_W <- boot_breakdown_W %>% filter(!is.na(var2))
	}

	for(s in unique(stratAssign_comp$stratum)){ # for each stratum
		for(sg in unique(fallback_rates[[1]]$stockGroup)){
			fStrat <- strataMatchUp %>% filter(stratum == s, stockGroup = sg) %>% pull(fallback)
			fBackRate <- 1 - (fallback_rates[[1]] %>% filter(stratum == fStrat, stockGroup == sg) %>% pull(p_fa))
			# selecting column that contains bootstraps for relevant fallback rate
			fBack_boot <- tibble(boot = (1:boots),
										fBack = 1 - fallback_rates[[2]][,which(fallback_rates[[1]]$stratum == fStrat &
																								fallback_rates[[1]]$stockGroup == sg)]
			)
			# H
			# groups that fall within this stockGroup
			tempGroups <- group_tables[[1]] %>% filter(stockGroup == sg) %>% select(-stockGroup)
			# tBool is boolean to select rows to adjust with the current value
			if(split_H_fallback == "var1"){
				tBool <- full_breakdown_H$stratum == s & full_breakdown_H$var1 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_H$stratum == s & boot_breakdown_H$var1 %in% unique(tempGroups[[1]])
			} else if (split_H_fallback == "var2"){
				tBool <- full_breakdown_H$stratum == s & full_breakdown_H$var2 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_H$stratum == s & boot_breakdown_H$var2 %in% unique(tempGroups[[1]])
			} else {
				rowIndex <- full_breakdown_H %>% mutate(index = 1:nrow(full_breakdown_H)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool <- (1:nrow(full_breakdown_H)) %in% rowIndex
				rowIndex <- boot_breakdown_H %>% mutate(index = 1:nrow(boot_breakdown_H)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool_boot <- (1:nrow(boot_breakdown_H)) %in% rowIndex
			}
			full_breakdown_H$total[tBool] <- full_breakdown_H$total[tBool] * fBackRate
			boot_breakdown_H$total[tBool_boot] <- boot_breakdown_H[tBool_boot,] %>%
				left_join(fBack_boot, by = "boot") %>% mutate(total = total * fBack) %>% pull(total)

			# HNC
			tempGroups <- group_tables[[2]] %>% filter(stockGroup == sg) %>% select(-stockGroup)
			if(split_HNC_fallback == "var1"){
				tBool <- full_breakdown_HNC$stratum == s & full_breakdown_HNC$var1 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_HNC$stratum == s & boot_breakdown_HNC$var1 %in% unique(tempGroups[[1]])
			} else if (split_HNC_fallback == "var2"){
				tBool <- full_breakdown_HNC$stratum == s & full_breakdown_HNC$var2 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_HNC$stratum == s & boot_breakdown_HNC$var2 %in% unique(tempGroups[[1]])
			} else {
				rowIndex <- full_breakdown_HNC %>% mutate(index = 1:nrow(full_breakdown_HNC)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool <- (1:nrow(full_breakdown_HNC)) %in% rowIndex
				rowIndex <- boot_breakdown_HNC %>% mutate(index = 1:nrow(boot_breakdown_HNC)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool_boot <- (1:nrow(boot_breakdown_HNC)) %in% rowIndex
			}
			full_breakdown_HNC$total[tBool] <- full_breakdown_HNC$total[tBool] * fBackRate
			boot_breakdown_HNC$total[tBool_boot] <- boot_breakdown_HNC[tBool_boot,] %>%
				left_join(fBack_boot, by = "boot") %>% mutate(total = total * fBack) %>% pull(total)

			# W
			tempGroups <- group_tables[[3]] %>% filter(stockGroup == sg) %>% select(-stockGroup)
			if(split_W_fallback == "var1"){
				tBool <- full_breakdown_W$stratum == s & full_breakdown_W$var1 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_W$stratum == s & boot_breakdown_W$var1 %in% unique(tempGroups[[1]])
			} else if (split_W_fallback == "var2"){
				tBool <- full_breakdown_W$stratum == s & full_breakdown_W$var2 %in% unique(tempGroups[[1]])
				tBool_boot <- boot_breakdown_W$stratum == s & boot_breakdown_W$var2 %in% unique(tempGroups[[1]])
			} else {
				rowIndex <- full_breakdown_W %>% mutate(index = 1:nrow(full_breakdown_W)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool <- (1:nrow(full_breakdown_W)) %in% rowIndex
				rowIndex <- boot_breakdown_W %>% mutate(index = 1:nrow(boot_breakdown_W)) %>% select(stratum, index, var1, var2) %>%
					filter(stratum == s) %>% inner_join(tempGroups, by = c("var1", "var2")) %>% pull(index)
				tBool_boot <- (1:nrow(boot_breakdown_W)) %in% rowIndex
			}
			full_breakdown_W$total[tBool] <- full_breakdown_W$total[tBool] * fBackRate
			boot_breakdown_W$total[tBool_boot] <- boot_breakdown_W[tBool_boot,] %>%
				left_join(fBack_boot, by = "boot") %>% mutate(total = total * fBack) %>% pull(total)
		}

	}

	# now calculate CIs for the run
	cis <- bind_rows(boot_breakdown_H, boot_breakdown_HNC, boot_breakdown_W) %>% group_by(rear, var1, var2, boot) %>%
		summarise(total = sum(total)) %>% ungroup %>% group_by(rear, var1, var2) %>%
		summarise(lci = quantile(total, alpha_ci/2), uci = quantile(total, 1 - (alpha_ci/2)))
	# calculate point estimates for the run and add CIs
	output <- bind_rows(full_breakdown_H, full_breakdown_HNC, full_breakdown_W) %>% group_by(rear, var1, var2) %>%
		summarise(pointEst = sum(total)) %>% full_join(cis, by = c("rear", "var1", "var2"))

	if(output_type == "summary") return(output = output)
	# mainly for troubleshooting or extending method
	if(output_type == "full") return(list(output = output, full_breakdown_H = full_breakdown_H,
													  full_breakdown_HNC = full_breakdown_HNC,
													  full_breakdown_W = full_breakdown_W,
													  boot_breakdown_H = boot_breakdown_H,
													  boot_breakdown_HNC = boot_breakdown_HNC,
													  boot_breakdown_W = boot_breakdown_W))

	# output W bootstraps by strata for pit array model
	W_boot <- bind_rows(boot_breakdown_H, boot_breakdown_HNC, boot_breakdown_W) %>% filter(rear == "W")
	# if var1 only estimates are included, prevent doubling the estimates
	if(any(is.na(W_boot$var2))) W_boot <- W_boot %>% filter(is.na(var2))
	W_boot <- W_boot %>% group_by(stratum, boot) %>% summarise(total = sum(total)) %>% spread(stratum, total)

	return(list(output = output, W_boot = W_boot))

}


# modular
# steps are: nighttime passage, fallback/reascension, composition of ascensions, multiply everything together
