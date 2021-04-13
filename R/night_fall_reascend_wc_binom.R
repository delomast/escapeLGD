# night time passage, fallback and reascension functions

#'
#' estimate fallback, reascension, and nighttime passage
#' @param full_reascend tibble with sWeek, stockGroup, numReascend, totalPass. sWeek is statistical week
#'   stockGroup is the stockGroup (fallback/reascension rates different between stock groups), numReascend is either the
#'   number of PITs ascending that were reascensions in that sWeek or the number of PITs ascending that later
#'   reascended (if using spillway data it MUST be the second)
#' @param full_night tibble with sWeek, nightPass, totalPass
#' @param stratAssign_fallback tibble with sWeek, stockGroup, stratum showing what stratum each sWeek
#'   corresponds to for each stockGroup
#' @param stratAssign_night tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for nighttime passage
#' @param boots The number of bootstap iterations to run
#' @param full_spillway tibble with "sWeek", "stockGroup", "laterAscend", "totalFall" columns. totalFall
#'   gives the number of PIT ascensions that were detected (with the spillway array) falling back and laterAscend
#'   gives the number of these that later reascended. If a strata has a summed value of 0 for totalFall,
#'   it will be assumed that all fallbacks result in reascension in that stratum. If a value of NULL is input
#'   for this, it will be assume that all fallbakcs result in reascension in all strata.
#'
#' @export
nightFall <- function(full_reascend, full_night, stratAssign_fallback, stratAssign_night,
							 boots = 2000, full_spillway = NULL){
	# some input checking
	if(ncol(full_reascend) != 4 || sum(c("sWeek", "stockGroup", "numReascend", "totalPass") %in%
		colnames(full_reascend)) != 4) stop("full_reascend must have 4 columns named: sWeek, stockGroup, numReascend, totalPass")
	if(ncol(full_night) != 3 || sum(c("sWeek", "nightPass", "totalPass") %in%
		colnames(full_night)) != 3) stop("full_night must have 3 columns named: sWeek, nightPass, totalPass")
	if(ncol(stratAssign_fallback) != 3 || sum(c("sWeek", "stockGroup", "stratum") %in%
		colnames(stratAssign_fallback)) != 3) stop("stratAssign_fallback must have 3 columns named: sWeek, stockGroup, stratum")
	if(ncol(stratAssign_night) != 2 || sum(c("sWeek", "stratum") %in%
		colnames(stratAssign_night)) != 2) stop("stratAssign_night must have 2 columns named: sWeek, stratum")
	if(!is.null(full_spillway) && (ncol(full_spillway) != 4 || sum(c("sWeek", "stockGroup", "laterAscend", "totalFall") %in%
		colnames(full_spillway)) != 2)) stop("full_spillway must have 4 columns named: sWeek, stockGroup, laterAscend, totalFall")


	# note that this allows different stratification for fallback by stockGroup, nighttime passage, and SCOBI
	# only constraint is that strata changes in fallback by stockGroup must coincide with strata changes in stratAssign_comp

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

	# make sure all totals are nonzero
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
			message("Assuming no fallback withOUT reascension for ", "stratum ",
					  fallback_rates[[1]]$stratum[i], " stockGroup ", fallback_rates[[1]]$stockGroup[i])
			fallback_rates[[1]]$p_fa[i] <- fallback_rates[[1]]$numReascend[i] / fallback_rates[[1]]$totalPass[i]
			fallback_rates[[2]][,i] <- rbinom(boots, fallback_rates[[1]]$totalPass[i],
														 fallback_rates[[1]]$p_fa[i]) / fallback_rates[[1]]$totalPass[i]
		} else { # otherwise infer both
			opts <- optim(par = c(.1,.9), fn = optimllh, gr = gradient_fallback_log_likelihood,
							  dfr = fallback_rates[[1]]$laterAscend[i], df = fallback_rates[[1]]$totalFall[i],
							  dr = fallback_rates[[1]]$numReascend[i], dt = fallback_rates[[1]]$totalPass[i],
							  control = list(fnscale = -1, maxit = 1000), method = "L-BFGS-B", upper = 1 - 1e-7,
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
#' @param nightPassage_rates the night passage rate part of the output of \code{nightFall}
#' @param wc tibble with two columns sWeek and wc. sWeek is statistical week and wc is the count of fish (NOT expanded for wc_prop)
#' @param wc_prop the proportion of the time fish are counted (e.g. 5/6)
#' @param stratAssign_night tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for nighttime passage
#' @param stratAssign_comp tibble with sWeek, stratum showing what stratum each sWeek corresponds to
#'   for composition estimation using the trap data
#' @param boots The number of bootstap iterations to run
#' @param alpha_ci The alpha value to use for calculating the CI of overall nighttime passage rate
#' @export
expand_wc_binom_night <- function(nightPassage_rates, wc, wc_prop, stratAssign_comp, stratAssign_night,
											 boots = 2000, alpha_ci = .1){
	# some input checking
	if(ncol(wc) != 2 || sum(c("sWeek", "wc") %in%
												  colnames(wc)) != 2) stop("wc must have 2 columns named: sWeek, wc")
	if(ncol(stratAssign_night) != 2 || sum(c("sWeek", "stratum") %in%
														colnames(stratAssign_night)) != 2) stop("stratAssign_night must have 2 columns named: sWeek, stratum")
	if(ncol(stratAssign_comp) != 2 || sum(c("sWeek", "stratum") %in%
														colnames(stratAssign_comp)) != 2) stop("stratAssign_comp must have 2 columns named: sWeek, stratum")
	# first, expand each sWeek for wc_prop and nighttime passage
	pointNight <- 0
	nightPass <- rep(0, boots)
	# similar setup for window count
	wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop))) # expanding for wc_prop
	wc_binom[[2]] <- matrix(nrow = boots, ncol = nrow(wc_binom[[1]]))
	for(i in 1:nrow(wc_binom[[1]])){
		# handle negative window counts
		if(wc_binom[[1]]$wc[i] < 0){
			wc_binom[[2]][,i] <- -rbinom(boots, -wc_binom[[1]]$wc[i], wc_prop) / wc_prop
		} else {
			wc_binom[[2]][,i] <- rbinom(boots, wc_binom[[1]]$wc[i], wc_prop) / wc_prop
		}
		# nighttime passage
		temp_nightStrata <- stratAssign_night$stratum[stratAssign_night$sWeek == wc$sWeek[i]]
		temp_nightStrata <- which(nightPassage_rates[[1]]$stratum == temp_nightStrata) # change into row/column number
		temp_nightRate <- nightPassage_rates[[1]]$p_night[temp_nightStrata]
		# expanding point estimate
		temp <- wc_binom[[1]]$wc[i] / (1 - temp_nightRate)
		pointNight <- pointNight + temp - wc_binom[[1]]$wc[i]
		wc_binom[[1]]$wc[i] <- temp
		rm(temp)
		# expanding bootstrap
		temp <- wc_binom[[2]][,i] / (1 - nightPassage_rates[[2]][,temp_nightStrata])
		nightPass <- nightPass + temp - wc_binom[[2]][,i]
		wc_binom[[2]][,i] <- temp
		rm(temp)
	}
	pointNight <- pointNight / sum(wc_binom[[1]]$wc)
	nightPass <- nightPass / rowSums(wc_binom[[2]])

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

	# add point estimate and CI for overall nighttime passage rate
	wc_binom[[3]] <- tibble(overallNightPass = pointNight,
									lci = quantile(nightPass, alpha_ci / 2),
									uci = quantile(nightPass, 1 - (alpha_ci / 2)))

	return(wc_binom)
}
