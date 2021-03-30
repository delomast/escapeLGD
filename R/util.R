# general utilities and other misc functions

#' A function to check if the composition and fallback strata are
#' compatible. Each composition stratum must have one and only one fallback stratum
#'
#' @param stratAssign_comp the strata assignment tibble for composition estimates
#' @param stratAssign_fallback the strata assignment tibble for fallback/reascension estimates
#' @param quiet when FALSE, it prints a message confirming that strata are compatible
#' @export
checkStrata <- function(stratAssign_comp, stratAssign_fallback, quiet = FALSE){
	eval <- stratAssign_comp %>% rename(comp = stratum) %>%
		full_join(stratAssign_fallback %>% rename(fall = stratum), by = "sWeek") %>%
		select(stockGroup, sWeek, comp, fall) %>% arrange(stockGroup, comp, fall)
	if(any(is.na(eval$comp))) stop("at least one sWeek in fallback strata that isn't present in composition strata")
	if(any(is.na(eval$fall))) stop("at least one sWeek in composition strata that isn't present with one or more stockGroup in fallback strata")
	for(g in unique(eval$stockGroup)){
		for(s in unique(eval$comp)){ # for each composition strata
			if((eval %>% filter(stockGroup == g, comp == s) %>% pull(fall) %>% n_distinct()) > 1) {
				stop("For stockGroup ", g, " in composition stratum ", s, " there is more than one fallback stratum")
			}
		}
	}
	if(!quiet) message("Strata are compatible")
}

#' a function to split PBT release groups estimated by var1
#' according to PIT tag detections. Used for WA "release groups"
#' that are split between groups released below and above LGD.
#' @param est_comp The output of either \code{ascension_composition} or
#'   \code{HNC_expand} where PBT release group was estimated as var1 for H and HNC.
#' @param splitByPITinput A tibble with columns releaseGroupPBT, stockGroup,
#'   releaseGroupPIT (optional), detectPIT, and tagRatePIT. Each line defines a
#'   PIT tag release group. releaseGroupPBT is the var1 (within H and HNC) that the PIT
#'   release group belongs to. stockGroup defines the stockGroup it belongs to (typically "upper" or "lower")
#'   detectPIT is the number of PIT tags detected from this group (actual number NOT expanded number)
#'   and tagRatePIT is the PIT tagging rate for this group.
#' @export
splitByPIT <- function(est_comp, splitByPITinput){
	uComb <- est_comp[[1]] %>% filter(rear %in% c("H", "HNC")) %>% select(var1) %>% distinct
	nBoot <- est_comp[[2]] %>% pull(boot) %>% n_distinct()
	for(g in uComb$var1[uComb$var1 %in% splitByPITinput$releaseGroupPBT]){
		# point estimate
		tempSplit <- splitByPITinput %>% filter(releaseGroupPBT == g) %>% mutate(point = detectPIT / tagRatePIT) %>%
			group_by(stockGroup) %>% summarise(r = sum(point), .groups = "drop") %>% mutate(r = r / sum(r), var1 = g)
		toBind <- est_comp[[1]] %>% filter(rear %in% c("H", "HNC"), var1 == g) %>% full_join(tempSplit, by = "var1") %>%
			mutate(var1 = paste0(var1, "_", stockGroup), total = total * r) %>% select(stratum, rear, var1, var2, total)
		est_comp[[1]] <- bind_rows(est_comp[[1]] %>% filter(!(rear %in% c("H", "HNC") & var1 == g)), toBind)
		# bootstrap estimate
		tempSplit <- splitByPITinput %>% filter(releaseGroupPBT == g) %>% mutate(point = detectPIT / tagRatePIT)
		toBind <- tibble()
		for(i in 1:nrow(tempSplit)){
			toBind <- toBind %>% bind_rows(
				tibble(boot = 1:nBoot, stockGroup = tempSplit$stockGroup[i],
						 bootPoint = rbinom(nBoot, round(tempSplit$point[i]), tempSplit$tagRatePIT[i]) / tempSplit$tagRatePIT[i])
			)
		}
		toBind <- toBind %>% group_by(boot, stockGroup) %>% summarise(r = sum(bootPoint), .groups = "drop_last") %>%
			mutate(r = r / sum(r), var1 = g) %>% ungroup()
		toBind <- est_comp[[2]] %>% filter(rear %in% c("H", "HNC"), var1 == g) %>%
			full_join(toBind, by = c("boot", "var1")) %>%
			mutate(var1 = paste0(var1, "_", stockGroup), total = total * r) %>% select(boot, stratum, rear, var1, var2, total)
		est_comp[[2]] <- bind_rows(est_comp[[2]] %>% filter(!(rear %in% c("H", "HNC") & var1 == g)), toBind)
	}

	return(est_comp)
}
