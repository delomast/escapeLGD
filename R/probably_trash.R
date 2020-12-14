#' one stratum, one dataset composition estimation
#' @param trap data for ONE stratum
# need to output number of fish by rear type, var1, var2
# ascension_composition <- function(trap, H_vars, HNC_vars, W_vars, HNC_W_adjust = TRUE, pbt_var = NULL, tagRates = NULL){
#
# 	pClip <- sum(trap$LGDMarkAD == "AD") / nrow(trap)
# 	tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI")
# 	pPhys <- sum(tempTrap$physTag) / nrow(tempTrap)
# 	# split no phystag into PBT-only HNC and unmarked, untagged
# 	tempTrap <- tempTrap %>% filter(LGDMarkAD == "AI", !physTag)
# 	if(nrow(tempTrap) == 0){
# 		pPBTonly <- 0
# 	} else {
# 		if(!any(!is.na(tempTrap$pbtAssign))){
# 			warning("all AI non-physTag samples in strata ", s, " were not assessed by PBT. Assuming all are Unassigned.")
# 			tempTrap$pbtAssign <- FALSE
# 		}
# 		tempTrap <- tempTrap %>% filter(!is.na(pbtAssign))
# 		if(!HNC_W_adjust){
# 			# no PBT expansion
# 			pbtOnly_rates[[1]] <- pbtOnly_rates[[1]] %>% bind_rows(
# 				stratum = s, pPBTonly = sum(tempTrap$pbtAssign) / nrow(tempTrap),
# 				totalTrap = nrow(tempTrap))
# 			pbtOnly_rates[[2]][,i] <- rbinom(boots, pbtOnly_rates[[1]]$totalTrap[i],
# 														pbtOnly_rates[[1]]$pPBTonly[i]) / pbtOnly_rates[[1]]$totalTrap[i]
# 		} else {
# 			# PBT tag rate expansion
# 			numHNC_pbtOnly <- tempTrap %>% filter(pbtAssign) %>% left_join(tagRates, by = pbt_var) %>%
# 				mutate(exp = 1 / tagRate) %>% pull(exp) %>% sum
# 			# truncate wild at 0
# 			if(numHNC_pbtOnly > nrow(tempTrap)) numHNC_pbtOnly <- nrow(tempTrap)
# 			pbtOnly_rates[[1]] <- pbtOnly_rates[[1]] %>% bind_rows(
# 				stratum = s, pPBTonly = numHNC_pbtOnly / nrow(tempTrap),
# 				totalTrap = nrow(tempTrap))
# 			# save seed so can regenerate these samples again during composition estimation
# 			# this preserves correlation in the estimates (groups with low tag rates more present
# 			# when HNC is higher and W is lower)
# 			# note that .Random.seed cannot be NULL here b/c rbinom was called earlier
# 			seed_for_HNC_W <- c(seed_for_HNC_W, .Random.seed)
# 			for(b in 1:boots){
# 				boot_numHNC_pbtOnly <- tempTrap %>% sample_n(nrow(.), replace = TRUE) %>%
# 					filter(pbtAssign) %>% left_join(tagRates, by = pbt_var) %>%
# 					mutate(exp = 1 / tagRate) %>% pull(exp) %>% sum
# 				if(boot_numHNC_pbtOnly > nrow(tempTrap)) boot_numHNC_pbtOnly <- nrow(tempTrap)
# 				pbtOnly_rates[[2]][b,i] <- boot_numHNC_pbtOnly / nrow(tempTrap)
# 			}
# 		}
# 	}
#
#
# }
