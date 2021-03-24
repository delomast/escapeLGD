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
	group_tables[[1]] <- group_tables[[1]] %>% select(-stockGroup) %>%
		left_join(H_groups, by = colnames(group_tables[[1]] %>% select(-stockGroup)))
	group_tables[[2]] <- group_tables[[2]] %>% select(-stockGroup) %>%
		left_join(HNC_groups, by = colnames(group_tables[[2]] %>% select(-stockGroup)))
	group_tables[[3]] <- group_tables[[3]] %>% select(-stockGroup) %>%
		left_join(W_groups, by = colnames(group_tables[[3]] %>% select(-stockGroup)))
	for(i in 1:3){
		if(any(is.na(group_tables[[i]]$stockGroup))){
			message("Missing an entry for ", c("H", "HNC", "W")[i], ". Returning a list of partially filled in templates.")
			return(group_tables)
		}
	}

	# matching up strata
	checkStrata(stratAssign_comp = stratAssign_comp, stratAssign_fallback = stratAssign_fallback, quiet = TRUE)
	strataMatchUp <- stratAssign_comp %>% full_join(stratAssign_fallback %>% rename(fallback = stratum), by = "sWeek") %>%
		select(-sWeek) %>% distinct

	# now multiply fallback rates as appropriate

	full_breakdown_H <- full_breakdown %>% filter(rear == "H")
	full_breakdown_HNC <- full_breakdown %>% filter(rear == "HNC")
	full_breakdown_W <- full_breakdown %>% filter(rear == "W")
	rm(full_breakdown) # save memory
	boot_breakdown_H <- boot_breakdown %>% filter(rear == "H")
	boot_breakdown_HNC <- boot_breakdown %>% filter(rear == "HNC")
	boot_breakdown_W <- boot_breakdown %>% filter(rear == "W")
	boots <- n_distinct(boot_breakdown$boot)
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
			fStrat <- strataMatchUp %>% filter(stratum == s, stockGroup == sg) %>% pull(fallback)
			if(length(fStrat) > 1) stop("Inconclusive matchup of strata for composition and fallback.")
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
		summarise(total = sum(total), .groups = "drop") %>% group_by(rear, var1, var2) %>%
		summarise(lci = quantile(total, alpha_ci/2), uci = quantile(total, 1 - (alpha_ci/2)), .groups = "drop")
	# calculate point estimates for the run and add CIs
	output <- bind_rows(full_breakdown_H, full_breakdown_HNC, full_breakdown_W) %>% group_by(rear, var1, var2) %>%
		summarise(pointEst = sum(total), .groups = "drop") %>% full_join(cis, by = c("rear", "var1", "var2"))

	# point estimates and CIs for rear type
	rearCIs <- bind_rows(boot_breakdown_H, boot_breakdown_HNC, boot_breakdown_W) %>% group_by(rear, var1, var2, boot) %>%
		summarise(total = sum(total), .groups = "drop") %>% group_by(rear, var1, var2) %>%
		summarise(lci = quantile(total, alpha_ci/2), uci = quantile(total, 1 - (alpha_ci/2)), .groups = "drop")
	rearCIs <- bind_rows(boot_breakdown_H, boot_breakdown_HNC, boot_breakdown_W)
	# avoid doubling estimates
	rearCIs <- tibble()
	rearPoint <- tibble()
	if(any(is.na(boot_breakdown_H$var2))) {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_H %>% filter(is.na(var2)))
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_H %>% filter(is.na(var2)))
	} else {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_H)
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_H)
	}
	if(any(is.na(boot_breakdown_HNC$var2))) {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_HNC %>% filter(is.na(var2)))
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_HNC %>% filter(is.na(var2)))

	} else {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_HNC)
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_HNC)

	}
	if(any(is.na(boot_breakdown_W$var2))) {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_W %>% filter(is.na(var2)))
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_W %>% filter(is.na(var2)))

	} else {
		rearCIs <- rearCIs %>% bind_rows(boot_breakdown_W)
		rearPoint <- rearPoint %>% bind_rows(full_breakdown_W)
	}
	rearCIs <- rearCIs  %>% group_by(rear, boot) %>%
		summarise(total = sum(total), .groups = "drop") %>% group_by(rear) %>%
		summarise(lci = quantile(total, alpha_ci/2), uci = quantile(total, 1 - (alpha_ci/2)), .groups = "drop")
	rear <- rearPoint %>% group_by(rear) %>% summarise(pointEst = sum(total), .groups = "drop") %>%
		full_join(rearCIs, by = "rear") %>% ungroup()

	if(output_type == "summary") return(list(output = output, rearType = rear))


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
	W_boot <- W_boot %>% group_by(stratum, boot) %>% summarise(total = sum(total), .groups = "drop") %>% spread(stratum, total)

	return(list(output = output, W_boot = W_boot))

}
