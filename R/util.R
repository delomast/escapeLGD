# general utilities and other misc functions

#' A function to check if the composition and fallback strata are
#' compatible. Each composition stratum must have one and only one fallback stratum
#'
#' @param stratAssign_comp the strata assignment tibble for composition estimates
#' @param stratAssign_fallback the strata assignment tibble for fallback/reascension estimates
#' @export
checkStrata <- function(stratAssign_comp, stratAssign_fallback){
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
	message("Strata are compatible")
}
