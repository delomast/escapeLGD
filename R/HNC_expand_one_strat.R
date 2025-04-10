# functions for calculating composition in one stratum expanding the HNC into the unmarked, untagged


#' one stratum, one dataset composition estimation, MLE
#' treating phystag and phys-untag separately
#' stepwise to deal with consistency and missing data issues
#' @param trap data for ONE stratum
#' @export
HNC_expand_one_strat_MLE <- function(trap, H_vars, HNC_vars, W_vars, wc_expanded,
											pbt_var = NULL, tagRates = NULL){
	pClip <- sum(trap$LGDMarkAD == "AD") / nrow(trap)
	tempTrap <- trap %>% filter(LGDMarkAD == "AI")
	pPhys <- sum(tempTrap$physTag) / nrow(tempTrap) # phys tag is FALSE if not assessed (NOT NA)

	# split no phystag into PBT-only HNC and unmarked, untagged
	tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI", !physTag)
	if(nrow(tempTrap) == 0){
		unmark_pbt_prop <- tibble()
		pPBTonly <- 0
	} else {
		if(!any(!is.na(tempTrap$pbtAssign))){
			warning("all AI non-physTag samples in strata ", s, " were not assessed by PBT. Assuming all are Unassigned.")
			tempTrap$pbtAssign <- FALSE
		}
		tempTrap <- tempTrap %>% filter(!is.na(pbtAssign))
		# PBT expansion
		#
		# will use estimates of pbtgroups vs wild comp in unmarked fish for both HNC and W composition
		# estimates. Estimating them here to make consistent across all without a joint estiamtion of everything
		# estimate pbt groups vs wild comp
		unmark_pbt_prop <- PBT_expand_calc_MLE(tempTrap[[pbt_var]], tagRates)
		pPBTonly <- unmark_pbt_prop %>% filter(group != "Unassigned") %>% pull(prop) %>% sum
	}

	estimates <- tibble()
	# H
	v1 <- H_vars[1]
	v2 <- if(length(H_vars) == 2) H_vars[2] else NULL
	# var1
	tempTrap <- trap %>% filter(LGDMarkAD == "AD", !is.na(trap[[v1]]))

	# need to skip iter or throw error if there are no H fish in tempTrap, but pClip is not 0
	if(nrow(tempTrap) == 0 && pClip > 0){
		# print("testing H NULL")
		return(NULL)
	}
	if(nrow(tempTrap) > 0){ # make sure there are H fish
		if(v1 == pbt_var){
			estim <- PBT_expand_calc_MLE(tempTrap[[v1]], tagRates)
		} else {
			props <- table(tempTrap[[v1]])
			estim <- tibble(group = names(props), prop = as.numeric(props/sum(props)))
		}
		# multiplying by expanded window count and other proportions then organizing output
		estim <- estim %>% mutate(rear = "H", var2 = NA, prop = prop * wc_expanded * pClip) %>%
			rename(var1 = group, total = prop) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)

		# var2
		if(!is.null(v2)){
			tempTrap <- tempTrap %>% filter(!is.na(tempTrap[[v2]]))
			tempEstim <- tibble()
			v2Data <- tibble(var1 = tempTrap[[v1]], var2 = tempTrap[[v2]]) %>%
				count(var1, var2)

			# make sure there is data
			sampleCount <- estim %>% filter(total > 0) %>% select(var1, total) %>% full_join(v2Data, by = "var1") %>%
				mutate(n = replace_na(n, 0)) %>% group_by(var1) %>%
				summarise(totNum = sum(n), .groups = "drop") %>% pull(totNum)
			# need to skip this iteration if bootstrapping, or throw an error if finding point estimate
			#  b/c no data to estimate composition of a subgroup (all samples are missing data)
			if(any(sampleCount < 1)){
				# print("testing H NULL v2")
				return(NULL)
			}

			if(v1 == pbt_var){
				if(!"Unassigned" %in% v2Data$var1 || !any(v2Data$var1 != "Unassigned")){
					#### handle when all fish assigned
					#### handle when all fish unassigned
					tempEstim <- v2Data %>% group_by(var1) %>% mutate(prop = n / sum(n), .groups = "drop")
				} else {
					#### handle when a mix of assigned and unassigned
					v2Data <- v2Data %>% complete(var1, var2, fill = list(n = 0)) %>%
						spread(var2, n)
					# put unassigned at the end
					v2Data <- bind_rows(v2Data %>% filter(var1 != "Unassigned"),
											  v2Data %>% filter(var1 == "Unassigned"))
					rn <- v2Data$var1
					v2Data <- v2Data %>% select(-var1) %>% as.matrix()
					rownames(v2Data) <- rn
					# v2Data is now matrix with counts, rows are var1, cols are var2
					piGroup <- tibble(group = rn) %>%
						left_join(estim %>% rename(group = var1) %>% mutate(prop = total/sum(total)), by = "group") %>%
						pull(prop)
					# piGroup is now a vector of "known" proportions of each group in order of rows
					tempTagRates <- tibble(group = rn) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

					# now run MLE routine with known piGroup
					varProbs <- PBT_var2_calc_MLE(v2Data = v2Data, piGroup = piGroup, tagRates = tempTagRates)
					for(i in 1:nrow(varProbs)){
						tempEstim <- bind_rows(tempEstim, tibble(var1 = rownames(varProbs)[i],
																			  var2 = colnames(varProbs),
																			  prop = varProbs[i,]))
					}
				}
			} else if (v2 == pbt_var){
				for(i in 1:nrow(estim)){
					v1Cat <- estim$var1[i]
					subGroupData <- tempTrap[[v2]][tempTrap[[v1]] == v1Cat]
					subEstim <- PBT_expand_calc_MLE(subGroupData, tagRates) %>%
						rename(var2 = group) %>% mutate(var1 = v1Cat) %>% select(var1, var2, prop)
					tempEstim <- bind_rows(tempEstim, subEstim)
				}
			} else {
				# neither is PBT var
				tempEstim <- v2Data %>% group_by(var1) %>% mutate(prop = n / sum(n), .groups = "drop")

			}

			estimates <- estimates %>% bind_rows(
				tempEstim %>% left_join(estim %>% select(var1, total), by = "var1") %>%
					mutate(total = prop * total, # multiply proportion by total number of fish
					rear = "H") %>% select(rear, var1, var2, total)
				)
		}
	}


	# HNC
	v1 <- HNC_vars[1]
	v2 <- if(length(HNC_vars) == 2) HNC_vars[2] else NULL
	# var1
	numHNC <- trap %>% filter(LGDMarkAD == "AI", (physTag | (!is.na(pbtAssign) & pbtAssign)),
										 !is.na(!!as.symbol(v1))) %>% nrow()
	tempTrap <- trap %>% filter(LGDMarkAD == "AI", !is.na(!!as.symbol(v1)))
	# need to skip iter or throw error if there are no HNC fish in tempTrap, but HNC estimate is not 0
	if(numHNC == 0 && (pPhys > 0 || !isTRUE(all.equal(pPBTonly, 0)))){
		# print("testing HNC NULL")
		return(NULL)
	}

	if(numHNC > 0){ # make sure there are HNC fish
		estim <- tibble()
		if(v1 == pbt_var){
			# phystag
			# expand phystag into unassigned phystag
			ptTrap <- tempTrap %>% filter(physTag)
			if(nrow(ptTrap) > 0){
				estim <- PBT_expand_calc_MLE(ptTrap[[v1]], tagRates) %>%
					mutate(total = prop * wc_expanded * (1 - pClip) * pPhys)
				estim_phystag_v2 <- estim %>% mutate(rear = NA, var2 = NA) %>%
					rename(var1 = group) %>% select(rear, var1, var2, total)
			} else if(pPhys > 0){
				# missing data issues, skip this bootstrap or throw an error
				# print("testing HNC v1.1 NULL")
				return(NULL)
			}
			rm(ptTrap)
			# non - phystag
			# expand into unmarked untagged fish (already calculated above)
			if(nrow(unmark_pbt_prop) > 0){
				estim <- bind_rows(estim, unmark_pbt_prop %>% filter(group != "Unassigned") %>% # Unassigned here is wild
						mutate(prop = prop / sum(prop),
								 total = prop * wc_expanded * (1 - pClip) * (1 - pPhys) * pPBTonly)
					) %>% group_by(group) %>% summarise(total = sum(total), .groups = "drop") %>% select(group, total)
			} else if(((1 - pPhys) * pPBTonly) > 0){
				# missing data issues, skip this bootstrap or throw an error
				return(NULL)
			}
		} else {
			# need to expand pbt only fish
			# don't need to expand physically tagged fish
			# phystag
			ptTrap <- tempTrap %>% filter(physTag)
			if(nrow(ptTrap) > 0){
				estim <- ptTrap %>% count(!!as.symbol(v1)) %>% mutate(prop = n / sum(n)) %>%
					rename(group = !!as.symbol(v1)) %>%
					mutate(total = prop * wc_expanded * (1 - pClip) * pPhys) %>% select(group, total)
			} else if(pPhys > 0){
				# missing data issues, skip this bootstrap or throw an error
				# print("testing HNC v1.2 NULL")
				return(NULL)
			}
			rm(ptTrap)
			# non - phystag
			# expand into unmarked untagged fish
			uptTrap <- tempTrap %>% filter(!physTag, !is.na(pbtAssign)) # so AI, no phystag, genotyped fish
			if(nrow(uptTrap) > 0 & any(uptTrap$pbtAssign)){
				v1Data <- tibble(var1 = uptTrap[[pbt_var]], var2 = uptTrap[[v1]]) %>%
					count(var1, var2)
				v1Data <- v1Data %>% complete(var1, var2, fill = list(n = 0)) %>%
					spread(var2, n)
				# put unassigned at the end
				v1Data <- bind_rows(v1Data %>% filter(var1 != "Unassigned"),
										  v1Data %>% filter(var1 == "Unassigned"))
				rn <- v1Data$var1
				v1Data <- v1Data %>% select(-var1) %>% as.matrix()
				rownames(v1Data) <- rn
				# v1Data is now matrix with counts, rows are var1, cols are var2
				piGroup <- tibble(group = rn) %>% left_join(unmark_pbt_prop, by = "group") %>%
					mutate(prop = prop / sum(prop)) %>% pull(prop)
				# piGroup is now a vector of "known" proportions of each group in order of rows
				tempTagRates <- tibble(group = rn) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

				# now run MLE routine with known piGroup
				varProbs <- PBT_var2_calc_MLE(v2Data = v1Data, piGroup = piGroup, tagRates = tempTagRates)
				tempEstim <- tibble()
				for(i in 1:nrow(varProbs)){
					tempEstim <- bind_rows(tempEstim, tibble(group = rownames(varProbs)[i],
																		  var1 = colnames(varProbs),
																		  prop_var2 = varProbs[i,]))
				}
				estim <- bind_rows(estim,
						tempEstim %>% left_join(unmark_pbt_prop, by = "group") %>%
							mutate(total = prop_var2 * prop * wc_expanded * (1 - pClip) * (1 - pPhys)) %>%
							filter(group != "Unassigned") %>% # Unassigned are wild here
							group_by(var1) %>% summarise(total = sum(total), .groups = "drop") %>% select(var1, total) %>%
							rename(group = var1)
					) %>% group_by(group) %>% summarise(total = sum(total), .groups = "drop")
			} else if(!isTRUE(all.equal(((1 - pPhys) * pPBTonly), 0))){
				# missing data issues, skip this bootstrap or throw an error
				# print("testing HNC v1.3 NULL")
				return(NULL)
			}
			rm(uptTrap)
		}
		estim <- estim %>% mutate(rear = "HNC", var2 = NA) %>%
			rename(var1 = group) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)


		# var2
		if(!is.null(v2)){
			if(v1 != pbt_var) stop("MLE is not currently an option for estimating the composition of HNC with two variables unless the first is the pbt_var")

			##### begin copy and paste


			# need to calculate comp for phystag and unphystag separately then add together
			# easiest to implement when PBT group is first variable since total and unphystag are already calculated

			# first check for data in each category
			tempTrap_check <- tempTrap %>% filter(physTag | pbtAssign, !is.na(!!as.symbol(v2)))
			tempv2Data <- tibble(var1 = tempTrap_check[[v1]], var2 = tempTrap_check[[v2]]) %>%
				count(var1, var2)
			rm(tempTrap_check)
			sampleCount <- estim %>% filter(total > 0) %>% select(var1, total) %>% full_join(tempv2Data, by = "var1") %>%
				mutate(n = replace_na(n, 0)) %>% group_by(var1) %>%
				summarise(totNum = sum(n), .groups = "drop") %>% pull(totNum)
			# need to skip this iteration if bootstrapping, or throw an error if finding point estimate
			#  b/c no data to estimate composition of a subgroup (all samples are missing data)
			if(any(sampleCount < 1)){
				# print("testing HNC NULL v2")
				return(NULL)
			}

			ptTrap <- tempTrap %>% filter(physTag, !is.na(!!as.symbol(v2)))
			tempEstim <- tibble()
			if(nrow(ptTrap) > 0){
				v2Data <- tibble(var1 = ptTrap[[v1]], var2 = ptTrap[[v2]]) %>%
					count(var1, var2)
				if(!"Unassigned" %in% v2Data$var1 || !any(v2Data$var1 != "Unassigned")){
					#### handle when all fish assigned
					#### handle when all fish unassigned
					tempEstim <- v2Data %>% group_by(var1) %>% mutate(prop = n / sum(n), .groups = "drop")
				} else {
					#### handle when a mix of assigned and unassigned
					v2Data <- v2Data %>% complete(var1, var2, fill = list(n = 0)) %>%
						spread(var2, n)
					# put unassigned at the end
					v2Data <- bind_rows(v2Data %>% filter(var1 != "Unassigned"),
											  v2Data %>% filter(var1 == "Unassigned"))
					rn <- v2Data$var1
					v2Data <- v2Data %>% select(-var1) %>% as.matrix()
					rownames(v2Data) <- rn
					# v2Data is now matrix with counts, rows are var1, cols are var2
					piGroup <- tibble(group = rn) %>%
						left_join(estim_phystag_v2 %>% rename(group = var1) %>%
									 	mutate(prop = total/sum(total)), by = "group") %>%
						pull(prop)
					# piGroup is now a vector of "known" proportions of each group in order of rows
					tempTagRates <- tibble(group = rn) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

					# now run MLE routine with known piGroup
					varProbs <- PBT_var2_calc_MLE(v2Data = v2Data, piGroup = piGroup, tagRates = tempTagRates)
					for(i in 1:nrow(varProbs)){
						tempEstim <- bind_rows(tempEstim, tibble(var1 = rownames(varProbs)[i],
																			  var2 = colnames(varProbs),
																			  prop = varProbs[i,]))
					}
				}

				tempEstim <- tempEstim %>% left_join(estim_phystag_v2 %>% select(var1, total), by = "var1") %>%
						mutate(total = prop * total # multiply proportion by total number of fish
							) %>% select(var1, var2, total)
			} else if(pPhys > 0){
				# missing data issues, skip this bootstrap or throw an error
				# print("testing HNC v1.2 NULL")
				return(NULL)
			}
			rm(ptTrap)
			# non - phystag
			# expand into unmarked untagged fish
			uptTrap <- tempTrap %>% filter(!physTag, !is.na(pbtAssign), !is.na(!!as.symbol(v2))) # so AI, no phystag, genotyped fish
			if(nrow(uptTrap) > 0 & any(uptTrap$pbtAssign)){
				v1Data <- tibble(var1 = uptTrap[[pbt_var]], var2 = uptTrap[[v2]]) %>%
					count(var1, var2)
				v1Data <- v1Data %>% complete(var1, var2, fill = list(n = 0)) %>%
					spread(var2, n)
				# put unassigned at the end
				v1Data <- bind_rows(v1Data %>% filter(var1 != "Unassigned"),
										  v1Data %>% filter(var1 == "Unassigned"))
				rn <- v1Data$var1
				v1Data <- v1Data %>% select(-var1) %>% as.matrix()
				rownames(v1Data) <- rn
				# v1Data is now matrix with counts, rows are var1, cols are var2
				piGroup <- tibble(group = rn) %>% left_join(unmark_pbt_prop, by = "group") %>%
					mutate(prop = prop / sum(prop)) %>% pull(prop)
				# piGroup is now a vector of "known" proportions of each group in order of rows
				tempTagRates <- tibble(group = rn) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

				# now run MLE routine with known piGroup
				varProbs <- PBT_var2_calc_MLE(v2Data = v1Data, piGroup = piGroup, tagRates = tempTagRates)
				tempEstim2 <- tibble()
				for(i in 1:nrow(varProbs)){
					tempEstim2 <- bind_rows(tempEstim2, tibble(group = rownames(varProbs)[i],
																		  var2 = colnames(varProbs),
																		  prop_var2 = varProbs[i,]))
				}
				tempEstim <- bind_rows(tempEstim,
										 tempEstim2 %>% left_join(unmark_pbt_prop, by = "group") %>%
										 	mutate(total = prop_var2 * prop * wc_expanded * (1 - pClip) * (1 - pPhys)) %>%
										 	filter(group != "Unassigned") %>% rename(var1 = group) %>%
										 	select(var1, var2, total)
				) %>% group_by(var1, var2) %>% summarise(total = sum(total), .groups = "drop")
			} else if(!isTRUE(all.equal(((1 - pPhys) * pPBTonly), 0))){
				# missing data issues, skip this bootstrap or throw an error
				# print("testing HNC v1.3 NULL")
				return(NULL)
			}
			rm(uptTrap)

			estimates <- estimates %>% bind_rows(
				tempEstim %>% mutate(rear = "HNC") %>% select(rear, var1, var2, total))
		}
	}


	# W
	v1 <- W_vars[1]
	v2 <- if(length(W_vars) == 2) W_vars[2] else NULL
	# var1
	tempTrap <- trap %>% filter(LGDMarkAD == "AI", !physTag, !is.na(pbtAssign),
										 !is.na(!!as.symbol(v1)))
	numW <- tempTrap %>% filter(!pbtAssign) %>% nrow()

	# need to skip iter or throw error if there are no W fish in tempTrap, but wild proportion is not 0
	if(numW == 0 && !isTRUE(all.equal((1 - pClip) * (1 - pPhys) * (1 - pPBTonly), 0))){
		# print("testing W NULL")
		return(NULL)
	}

	if(numW > 0){ # make sure there are W fish
		estim <- tibble()
		v1Data <- tibble(var1 = tempTrap[[pbt_var]], var2 = tempTrap[[v1]]) %>%
			count(var1, var2)

		if(!any(v1Data$var1 != "Unassigned")){
			##### handle if all Unassigned
			tempEstim <- v1Data %>% select(var2, n) %>% mutate(prop = n / sum(n)) %>%
				select(var2, prop) %>% rename(var1 = var2) # renaming b/c pbt variable is removed
		} else {
			#### handle if mix
			v1Data <- v1Data %>% complete(var1, var2, fill = list(n = 0)) %>%
				spread(var2, n)
			# put unassigned at the end
			v1Data <- bind_rows(v1Data %>% filter(var1 != "Unassigned"),
									  v1Data %>% filter(var1 == "Unassigned"))
			rn <- v1Data$var1
			v1Data <- v1Data %>% select(-var1) %>% as.matrix()
			rownames(v1Data) <- rn
			# v1Data is now matrix with counts, rows are var1, cols are var2
			piGroup <- tibble(group = rn) %>% left_join(unmark_pbt_prop, by = "group") %>%
				mutate(prop = prop / sum(prop)) %>% pull(prop)
			# piGroup is now a vector of "known" proportions of each group in order of rows
			tempTagRates <- tibble(group = rn) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

			# now run MLE routine with known piGroup
			varProbs <- PBT_var2_calc_MLE(v2Data = v1Data, piGroup = piGroup, tagRates = tempTagRates)
			tempEstim <- tibble(var1 = colnames(varProbs),
									  prop = varProbs["Unassigned",]) %>%
				# MLE can't push parameters all the way to zero, and non-zero estimates for unsampled groups
				#  can cause "NULL" to be returned during estimates for var2. So manaully settign any group
				#  not present in the untagged - unassigned fish to 0
				left_join(tibble(var1 = colnames(v1Data), numTrap = v1Data["Unassigned",]), by = "var1") %>%
				filter(numTrap > 0) %>% mutate(prop = prop / sum(prop)) %>% select(var1, prop)
		}

		estim <- tempEstim %>% mutate(total = prop * wc_expanded * (1 - pClip) * (1 - pPhys) * (1 - pPBTonly)) %>%
								 	mutate(rear = "W", var2 = NA) %>% select(rear, var1, var2, total)
		estimates <- estimates %>% bind_rows(estim)


		# var2
		if(!is.null(v2)){
			###
			# this is the tricky one. have to fix group proportions and var1 proportions
			# then infer proportions of var2 within var1
			# going to assume var1 and var2 are independent in HNC? Yes b/c var1 should be genetic stock
			###

			estim2 <- tibble()
			tempTrap <- tempTrap %>% filter(!is.na(!!as.symbol(v2)))

			# make sure there is data
			sampCounts <- estim %>% select(var1) %>% left_join(tempTrap %>%
				filter(!pbtAssign) %>% count(!!as.symbol(v1)) %>% rename(var1 = !!as.symbol(v1)), by = "var1") %>%
				mutate(n = replace_na(n, 0)) %>% pull(n)
			if(any(sampCounts < 1)){
				# print("testing W v2 NULL")
				return(NULL)
			}

			v2Data <- tibble(var1 = tempTrap[[pbt_var]], var2 = tempTrap[[v1]], var3 = tempTrap[[v2]]) %>%
				count(var1, var2, var3)

			if(!any(v2Data$var1 != "Unassigned")){
				##### handle if all Unassigned
				estim2 <- v2Data %>% select(var2, var3, n) %>% group_by(var2) %>% mutate(prop = n / sum(n)) %>%
					select(var2, var3, prop) %>% rename(var1 = var2, var2 = var3) # renaming b/c pbt variable is removed
			} else {
				#### handle if mix of assigned and unassigned
				v2Data <- v2Data %>% complete(var1, var2, var3, fill = list(n = 0))

				# now split into two for pbt assigned, one for Unassigned
				hnc_var3 <- v2Data %>% filter(var1 != "Unassigned") %>% group_by(var1, var3) %>%
					summarise(n = sum(n), .groups = "drop") %>% spread(var3, n)
				un_counts <- v2Data %>% filter(var1 == "Unassigned") %>% select(var2, var3, n) %>%
					spread(var3, n)
				# make into matrices and make orders the same
				rn_pbt <- hnc_var3$var1
				rn_v2 <- un_counts$var2

				un_counts <- un_counts %>% select(-var2) %>% as.matrix()
				rownames(un_counts) <- rn_v2
				hnc_var3 <- hnc_var3 %>% select(-var1) %>% as.matrix()
				rownames(hnc_var3) <- rn_pbt
				hnc_var3 <- hnc_var3[, colnames(un_counts), drop = FALSE] # make column orders the same

				piGroup <- tibble(group = c(rn_pbt, "Unassigned")) %>% left_join(unmark_pbt_prop, by = "group") %>%
					mutate(prop = prop / sum(prop)) %>% pull(prop)
				# piGroup is now a vector of "known" proportions of each group in order of rows
				tempTagRates <- tibble(group = rn_pbt) %>% left_join(tagRates, by = "group") %>% pull(tagRate)

				# rows are pbt groups + unassigned, columns are var2(of 3) categories in order of rows of un_counts
				varProbs <- varProbs[c(rn_pbt, "Unassigned"), rn_v2]

				# inputs:
				# knowns:
				# 	piGroup
				# 	tempTagRates
				# 	varProbs - known composition of var2(of 3)
				# data:
				# hnc_var3 - counts of PBT assigned in cats of var3
				# un_counts - counts of unassigned in cats of var2 (rows) and var3 (cols)
				var3Probs <- PBT_var3_calc_MLE(piGroup = piGroup, tagRates = tempTagRates,
														 var2Probs = varProbs, hnc_var3 = hnc_var3,
														 un_counts = un_counts)
				for(i in 1:nrow(var3Probs)){
					estim2 <- estim2 %>% bind_rows(tibble(var1 = rownames(var3Probs)[i],
																			  var2 = colnames(var3Probs),
																			  prop = var3Probs[i,]))
				}


			}

			estimates <- estimates %>% bind_rows(
				estim2 %>% left_join(estim %>% select(var1, total), by = "var1") %>%
					filter(!is.na(total)) %>% # some groups have composition estimated (stays at initial values) when present in HNC but not in W, need to remove b/c total is 0
					mutate(total = prop * total, # multiply proportion by total number of fish
							 rear = "W") %>% select(rear, var1, var2, total)
			)
		}
	}
	return(list(estimates = estimates,
					testing = list(pClip = pClip, pPhys = pPhys, pPBTonly = pPBTonly)))


}














###################################################################################################



#' one stratum, one dataset composition estimation, accounting
#' @param trap data for ONE stratum
#' @export
HNC_expand_one_strat <- function(trap, H_vars, HNC_vars, W_vars, wc_expanded,
											pbt_var = NULL, tagRates = NULL){

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
			props <- table(tempTrap[[v1]]) # removes NAs automatically
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
				tempEstimSwitch <- TRUE
				if(length(subGroupData) == 0){
					if(v1 == pbt_var && v1Cat == "Unassigned" && estim$total[i] == 0){
						# this can happen when all fish assign
						# let's just report all categories as 0 so user knows it wasn't left out
						tempEstim <- tibble(var2 = unique(tempTrap[[v2]]), prop = 0)
						tempEstimSwitch <- FALSE
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
					if(tempEstimSwitch){
						props <- table(subGroupData)
						tempEstim <- tibble(var2 = names(props), prop = as.numeric(props))
					}
					if(nrow(assignedData) > 0){
						tempEstim <- tempEstim %>% left_join(assignedData, by = "var2") %>%
							mutate(expectUnass = replace_na(expectUnass, 0), prop = prop - expectUnass)
						tempEstim$prop[tempEstim$prop < 0] <- 0 # truncate at 0
					}
					tempEstim <- tempEstim %>% mutate(prop = prop/sum(prop)) %>% rename(group = var2) %>% select(group, prop)
				} else {
					props <- table(subGroupData)
					tempEstim <- tibble(group = names(props), prop = as.numeric(props/sum(props)))
				}
				if(estim$total[i] == 0){
					tempEstim <- tempEstim %>% mutate(prop = 0, rear = "H", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				} else {
					tempEstim <- tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
																 rear = "H", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				}
				estim2 <- estim2 %>% bind_rows(tempEstim)
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
			estim$prop[estim$group == "Unassigned"] <- max(estim$prop[estim$group == "Unassigned"] - removeUnass, 0)
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
					tempEstim$prop[tempEstim$group == "Unassigned"] <- max(tempEstim$prop[tempEstim$group == "Unassigned"] - removeUnass, 0)
					tempEstim <- tempEstim %>% group_by(group) %>% summarise(prop = sum(prop), .groups = "drop") %>%
						mutate(prop = prop / sum(prop))
				} else if(v1 == pbt_var && v1Cat == "Unassigned"){
					# need to correct composition of Unassigned subGroup estimates for v2
					# first sample counts in Unassigned
					props <- table(subGroupData[[v2]])
					tempEstim <- tibble(group = names(props), prop = as.numeric(props))
					# now expected number of phystag unassigned due to tag rates
					tempSubtract_data <- tempTrap %>% filter(physTag, pbtAssign)
					if(nrow(tempSubtract_data) > 0){
						tempEstim <- tempEstim %>% left_join(tempSubtract_data %>%
							group_by(!!as.symbol(v2)) %>% summarise(diff = sum(diff), .groups = "drop") %>%
							rename(group = !!as.symbol(v2)) %>% select(group, diff), by = "group") %>%
							mutate(diff = replace_na(diff, 0), prop = prop - diff)
						tempEstim$prop[tempEstim$prop < 0] <- 0 # truncate at 0
					}
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
				if(estim$total[i] == 0){
					tempEstim <- tempEstim %>% mutate(prop = 0, rear = "HNC", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				} else {
					tempEstim <- tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
																 rear = "HNC", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				}
				estim2 <- estim2 %>% bind_rows(tempEstim)
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
		estim <- tibble(group = names(props), prop = as.numeric(props))
		if(nrow(tempTrap_HNC) > 0){
			estim <- estim %>% left_join(
				tempTrap_HNC %>% group_by(!!as.symbol(v1)) %>% summarise(diff = sum(diff), .groups = "drop") %>%
					rename(group = !!as.symbol(v1)),
				by = "group") %>% mutate(diff = replace_na(diff, 0), prop = prop - diff) # subtract expanded HNC
			estim$prop[estim$prop < 0] <- 0 # truncate at 0
		}
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
					props <- table(subGroupData[[v2]])
					tempEstim <- tibble(group = names(props), prop = as.numeric(props))
					if(nrow(tempTrap_HNC) > 0){
						tempEstim <- tempEstim %>% left_join(
							tempTrap_HNC %>% filter(!!as.symbol(v1) == v1Cat) %>% group_by(!!as.symbol(v2)) %>%
								summarise(diff = sum(diff), .groups = "drop") %>% rename(group = !!as.symbol(v2)),
							by = "group") %>% mutate(diff = replace_na(diff, 0), prop = prop - diff) # subtract expanded HNC
						tempEstim$prop[tempEstim$prop < 0] <- 0 # truncate at 0
					}
					tempEstim <- tempEstim %>% mutate(prop = prop / sum(prop))
				}

				if(estim$total[i] == 0){
					tempEstim <- tempEstim %>% mutate(prop = 0, rear = "W", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				} else {
					tempEstim <- tempEstim %>% mutate(prop = prop * estim$total[i], # multiply proportion by total number of fish
																 rear = "W", var1 = v1Cat) %>%
						rename(var2 = group, total = prop) %>% select(rear, var1, var2, total)
				}

				estim2 <- estim2 %>% bind_rows(tempEstim)
			}
			estimates <- estimates %>% bind_rows(estim2)
		}
	}
	return(list(estimates = estimates,
					testing = list(pClip = pClip, pPhys = pPhys, pPBTonly = pPBTonly)))
}
