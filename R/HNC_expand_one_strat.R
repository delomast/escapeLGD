# need to output number of fish by rear type, var1, var2

#' one stratum, one dataset composition estimation
#' @param trap data for ONE stratum
#' @export
HNC_expand_one_strat <- function(trap, H_vars, HNC_vars, W_vars, wc_expanded,
											pbt_var = NULL, tagRates = NULL){


	colnames(tagRates) <- c("group", "tagRate")
	if("Unassigned" %in% tagRates[[1]]){
		if(tagRates[[2]][tagRates[[1]] == "Unassigned"] != 1){
			stop("Unassigned must not be included in the tag rate file or have a tag rate of 1")
		}
	} else {
		# add Unassigned with tag rate of 1
		tagRates <- bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1))
	}

	pClip <- sum(trap$LGDMarkAD == "AD") / nrow(trap)
	tempTrap <- trap %>% filter(LGDMarkAD == "AI")
	pPhys <- sum(tempTrap$physTag) / nrow(tempTrap) # phys tag is FALSE if not assessed (NOT NA)

	# split no phystag into PBT-only HNC and unmarked, untagged
	tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI", !physTag)
	if(nrow(tempTrap) == 0){
		pPBTonly <- 0
	} else {
		if(!any(!is.na(tempTrap$pbtAssign))){
			warning("all AI non-physTag samples in strata ", s, " were not assessed by PBT. Assuming all are Unassigned.")
			tempTrap$pbtAssign <- FALSE
		}
		tempTrap <- tempTrap %>% filter(!is.na(pbtAssign))
		# PBT expansion
		pPBTonly <- PBT_expand_calc(tempTrap[[pbt_var]], tagRates) %>%
			filter(group != "Unassigned") %>% pull(prop) %>% sum
	}

	estimates <- tibble()
	# H
	v1 <- H_vars[1]
	v2 <- if(length(H_vars) == 2) H_vars[2] else NULL
	# var1
	tempTrap <- trap %>% filter(LGDMarkAD == "AD", !is.na(trap[[v1]]))

	# need to skip iter or throw error if there are no H fish in tempTrap, but pClip is not 0
	if(nrow(tempTrap) == 0 && pClip > 0) return(NULL)

	if(nrow(tempTrap) > 0){ # make sure there are H fish
		if(v1 == pbt_var){
			estim <- PBT_expand_calc(tempTrap[[v1]], tagRates)
		} else {
			props <- table(tempTrap[[v2]]) # removes NAs automatically
			estim <- tibble(group = names(props), prop = as.numeric(props/sum(props)))
		}
		# multiplying by expanded window count and other proportions then organizing output
		estim <- estim %>% mutate(rear = "H", var2 = NA, prop = prop * wc_expanded * pClip) %>%
			rename(var1 = group, total = prop) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)

		# var2
		if(!is.null(v2)){
			estim2 <- tibble()
			tempTrap <- tempTrap %>% filter(!is.na(tempTrap[[v2]]))
			# Unassigned may be 0 if pbt_var is var1
			for(i in 1:nrow(estim)){
				v1Cat <- estim$var1[i]
				subGroupData <- tempTrap[[v2]][tempTrap[[v1]] == v1Cat]
				if(length(subGroupData) == 0){
					if(v1 == pbt_var && v1Cat == "Unassigned" && estim$total[i] == 0){
						# this can happen when all fish assign
						# let's just report all categories as 0 so user knows it wasn't left out
						tempEstim <- tibble(group = unique(tempTrap[[v2]]), prop = 0)
					} else {
						# need to skip this iteration if bootstrapping, or throw an error if finding point estimate
						#  b/c no data to estimate composition of this subgroup (all samples are missing data)
						return(NULL)
					}
				}
				if(v2 == pbt_var){
					tempEstim <- PBT_expand_calc(subGroupData, tagRates)
				} else if(v1 == pbt_var && v1Cat == "Unassigned"){
					# need to correct composition of Unassigned subGroup estimates for v2
					assignedData <- tibble(group = tempTrap[[v1]][tempTrap[[v1]] != "Unassigned"],
												  var2 = tempTrap[[v2]][tempTrap[[v1]] != "Unassigned"]) %>%
						left_join(tagRates, by = "group") %>% mutate(expectUnass = (1 / tagRate) - 1) %>%
						group_by(var2) %>% summarise(expectUnass = sum(expectUnass), .groups = "drop") %>% select(var2, expectUnass)
					props <- table(subGroupData)
					tempEstim <- tibble(var2 = names(props), prop = as.numeric(props)) %>%
						left_join(assignedData, by = "var2") %>% mutate(expectUnass = replace_na(expectUnass, 0),
						prop = prop - expectUnass)
					tempEstim$prop[tempEstim$prop < 0] <- 0 # truncate at 0
					tempEstim <- tempEstim %>% mutate(prop = prop/sum(prop)) %>% rename(group = var2) %>% select(group, prop)
				} else {
					props <- table(subGroupData)
					tempEstim <- tibble(group = names(props), prop = as.numeric(props/sum(props)))
				}
				estim2 <- estim2 %>% bind_rows(
					tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
												rear = "H", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total))
			}
			estimates <- estimates %>% bind_rows(estim2)
		}
	}

	# HNC
	v1 <- HNC_vars[1]
	v2 <- if(length(HNC_vars) == 2) HNC_vars[2] else NULL
	# var1
	tempTrap <- trap %>% filter(LGDMarkAD == "AI", (physTag | (!is.na(pbtAssign) & pbtAssign)),
										 !is.na(!!as.symbol(v1)))
	tempTrap_HNC <- tibble() # save for use in W decomposition later, initialize in case no HNC fish

	# need to skip iter or throw error if there are no HNC fish in tempTrap, but HNC estimate is not 0
	if(nrow(tempTrap) == 0 && (pPhys > 0 || pPBTonly > 0)) return(NULL)

	if(nrow(tempTrap) > 0){ # make sure there are HNC fish
		# expanding and truncating as needed from the beginning
		tempTag <- tagRates
		colnames(tempTag)[1] <- pbt_var
		tempTrap <- tempTrap %>% left_join(tempTag, by = pbt_var) %>%
			mutate(tagRate = replace_na(tagRate, 1), exp = 1 / tagRate, diff = exp - 1)
		rm(tempTag)
		# truncating at number of unassigned fish for phystag fish
		max_expand <- sum(tempTrap$physTag & !is.na(tempTrap[[pbt_var]]) & tempTrap[[pbt_var]] == "Unassigned")
		if(sum(tempTrap$diff[tempTrap$physTag]) > max_expand){
			tempTrap$diff[tempTrap$physTag] <- max_expand * (tempTrap$diff[tempTrap$physTag] / sum(tempTrap$diff[tempTrap$physTag]))
		}
		# truncating at number of unmarked untagged fish for pbtonly fish
		max_expand <- sum(trap$LGDMarkAD == "AI" & !trap$physTag & !is.na(trap$pbtAssign) & !trap$pbtAssign)
		if(sum(tempTrap$diff[!tempTrap$physTag]) > max_expand){
			tempTrap$diff[!tempTrap$physTag] <- max_expand * (tempTrap$diff[!tempTrap$physTag] / sum(tempTrap$diff[!tempTrap$physTag]))
		}
		tempTrap$exp <- 1 + tempTrap$diff
		tempTrap_HNC <- tempTrap %>% filter(!physTag) # want to save PBT only fish for W decomposition

		estim <- tibble()
		if(v1 == pbt_var){
			estim <- tempTrap %>% group_by(!!as.symbol(pbt_var), physTag) %>%
				summarise(prop = sum(exp), diff = sum(diff), .groups = "drop") %>%
				rename(group = !!as.symbol(pbt_var)) %>% select(group, physTag, prop, diff)
			removeUnass <- estim %>% filter(physTag) %>% pull(diff) %>% sum
			if(!any(estim$group == "Unassigned")) estim <- estim %>% bind_rows(tibble(group = "Unassigned", physTag = TRUE, prop = 0, diff = 0))
			# expand phystag into unassigned phystag
			estim$prop[estim$group == "Unassigned"] <- estim$prop[estim$group == "Unassigned"] - removeUnass
			estim <- estim %>% group_by(group) %>% summarise(prop = sum(prop), .groups = "drop") %>% mutate(prop = prop / sum(prop))
		} else {
			# need to expand pbt only fish, truncating for number of wild fish
			# don't need to expand physically tagged fish
			estim <- tempTrap %>% group_by(!!as.symbol(v1), physTag) %>%
				summarise(prop = sum(exp), diff = sum(diff), .groups = "drop") %>%
				 rename(group = !!as.symbol(v1)) %>% select(group, physTag, prop, diff)
			estim$prop[estim$physTag] <- estim$prop[estim$physTag] - estim$diff[estim$physTag] # unexpanding phystag fish
			estim <- estim %>% group_by(group) %>% summarise(prop = sum(prop), .groups = "drop") %>% mutate(prop = prop / sum(prop))
		}
		# multiplying by expanded window count and other proportions then organizing output
		estim <- estim %>% mutate(rear = "HNC", var2 = NA, prop = prop * wc_expanded * (1 - pClip) *
										  	(pPhys + ((1 - pPhys) * pPBTonly))) %>%
			rename(var1 = group, total = prop) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)


		# var2
		if(!is.null(v2)){
			estim2 <- tibble()
			tempTrap <- tempTrap %>% filter(!is.na(!!as.symbol(v2)))
			# Unassigned may be 0 if pbt_var is var1
			for(i in 1:nrow(estim)){
				v1Cat <- estim$var1[i]
				subGroupData <- tempTrap %>% filter(!!as.symbol(v1) == v1Cat)
				if(nrow(subGroupData) == 0){
					if(v1 == pbt_var && v1Cat == "Unassigned" && estim$total[i] == 0){
						# this can happen when all fish assign
						# let's just report all categories as 0 so user knows it wasn't left out
						tempEstim <- tibble(group = unique(tempTrap[[v2]]), prop = 0)
					} else {
						# need to skip this iteration if bootstrapping, or throw an error if finding point estimate
						#  b/c no data to estimate composition of this subgroup (all samples are missing data)
						return(NULL)
					}
				} else if(v2 == pbt_var){
					# need to expand phystag and non separately
					tempEstim <- subGroupData %>% group_by(!!as.symbol(pbt_var), physTag) %>%
						summarise(prop = sum(exp), diff = sum(diff), .groups = "drop") %>%
						rename(group = !!as.symbol(pbt_var)) %>% select(group, physTag, prop, diff)
					removeUnass <- tempEstim %>% filter(physTag) %>% pull(diff) %>% sum
					if(!any(tempEstim$group == "Unassigned")) tempEstim <- tempEstim %>%
						bind_rows(tibble(group = "Unassigned", physTag = TRUE, prop = 0, diff = 0))
					# expand phystag into unassigned phystag
					tempEstim$prop[tempEstim$group == "Unassigned"] <- tempEstim$prop[tempEstim$group == "Unassigned"] - removeUnass
					tempEstim <- tempEstim %>% group_by(group) %>% summarise(prop = sum(prop), .groups = "drop") %>%
						mutate(prop = prop / sum(prop))
				} else if(v1 == pbt_var && v1Cat == "Unassigned"){
					# need to correct composition of Unassigned subGroup estimates for v2
					# first sample counts in Unassigned
					props <- table(subGroupData[[v2]])
					tempEstim <- tibble(group = names(props), prop = as.numeric(props))
					# now expected number of phystag unassigned due to tag rates
					tempEstim <- tempEstim %>% left_join(tempTrap %>% filter(physTag, pbtAssign) %>%
						group_by(!!as.symbol(v2)) %>% summarise(diff = sum(diff), .groups = "drop") %>%
							rename(group = !!as.symbol(v2)) %>% select(group, diff), by = "group") %>%
						mutate(diff = replace_na(diff, 0), prop = prop - diff)
					tempEstim$prop[tempEstim$prop < 0] <- 0 # truncate at 0
					tempEstim <- tempEstim %>% mutate(prop = prop / sum(prop)) %>% select(group, prop)
				} else if (v1 == pbt_var) {
					# don't need to expand anything b/c tagged vs untagged is random
					props <- table(subGroupData[[v2]])
					tempEstim <- tibble(group = names(props), prop = as.numeric(props/sum(props)))
				} else {
					# neither variable is PBT variable
					# need to expand pbt only fish, don't need to expand physically tagged fish
					tempEstim <- subGroupData %>% group_by(!!as.symbol(v2), physTag) %>%
						summarise(prop = sum(exp), diff = sum(diff), .groups = "drop") %>%
						rename(group = !!as.symbol(v2)) %>% select(group, physTag, prop, diff)
					tempEstim$prop[tempEstim$physTag] <- tempEstim$prop[tempEstim$physTag] -
						tempEstim$diff[tempEstim$physTag] # unexpanding phystag fish
					tempEstim <- tempEstim %>% group_by(group) %>% summarise(prop = sum(prop), .groups = "drop") %>%
						mutate(prop = prop / sum(prop))
				}

				estim2 <- estim2 %>% bind_rows(
					tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
												rear = "HNC", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total))
			}
			estimates <- estimates %>% bind_rows(estim2)
		}
	}


	# W
	v1 <- W_vars[1]
	v2 <- if(length(W_vars) == 2) W_vars[2] else NULL
	# var1
	tempTrap <- trap %>% filter(LGDMarkAD == "AI", !physTag, !is.na(pbtAssign), !pbtAssign,
										 !is.na(!!as.symbol(v1)))

	# need to skip iter or throw error if there are no W fish in tempTrap, but wild proportion is not 0
	if(nrow(tempTrap) == 0 && (1 - pClip) * (1 - pPhys) * (1 - pPBTonly) > 0) return(NULL)

	if(nrow(tempTrap) > 0){ # make sure there are W fish
		# tempTrap_HNC has PBT only fish with a value for v1
		# no PBT for W fish!
		props <- table(tempTrap[[v1]]) # removes NAs automatically
		estim <- tibble(group = names(props), prop = as.numeric(props)) %>% left_join(
			tempTrap_HNC %>% group_by(!!as.symbol(v1)) %>% summarise(diff = sum(diff), .groups = "drop") %>%
				rename(group = !!as.symbol(v1)),
			by = "group") %>% mutate(prop = prop - diff) # subtract expanded HNC
		estim$diff[estim$diff < 0] <- 0 # truncate at 0
		estim <- estim %>% mutate(prop = prop / sum(prop))

		# multiplying by expanded window count and other proportions then organizing output
		estim <- estim %>% mutate(rear = "W", var2 = NA, prop = prop * wc_expanded * (1 - pClip) *
										  	(1 - pPhys) * (1 - pPBTonly)) %>%
			rename(var1 = group, total = prop) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)

		# var2
		if(!is.null(v2)){
			estim2 <- tibble()
			tempTrap <- tempTrap %>% filter(!is.na(!!as.symbol(v2)))
			for(i in 1:nrow(estim)){
				v1Cat <- estim$var1[i]
				subGroupData <- tempTrap %>% filter(!!as.symbol(v1) == v1Cat)
				if(nrow(subGroupData) == 0){
					# need to skip this iteration if bootstrapping, or throw an error if finding point estimate
					#  b/c no data to estimate composition of this subgroup (all samples are missing data)
					return(NULL)
				} else {
					props <- table(tempTrap[[v2]])
					tempEstim <- tibble(group = names(props), prop = as.numeric(props)) %>% left_join(
						tempTrap_HNC %>% filter(!!as.symbol(v1) == v1Cat) %>% group_by(!!as.symbol(v2))
							%>% summarise(diff = sum(diff), .groups = "drop") %>% rename(group = !!as.symbol(v2)),
						by = "group") %>% mutate(prop = prop - diff) # subtract expanded HNC
					tempEstim$diff[tempEstim$diff < 0] <- 0 # truncate at 0
					tempEstim <- tempEstim %>% mutate(prop = prop / sum(prop))
				}

				estim2 <- estim2 %>% bind_rows(
					tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
												rear = "W", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total))
			}
			estimates <- estimates %>% bind_rows(estim2)
		}
	}
	return(list(estimates = estimates,
					testing = list(pClip = pClip, pPhys = pPhys, pPBTonly = pPBTonly)))
}
