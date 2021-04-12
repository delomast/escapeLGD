library(tidyverse)
library(lubridate)

numBoots <- 20
spp <- "STHD"
nf <- nightFall(full_reascend = fullReComplete, full_night = fullNi,
					 stratAssign_fallback = stratFallback, stratAssign_night = stratNight,
					 full_spillway = NULL, boots = numBoots)

exp_wc <- expand_wc_binom_night(nightPassage_rates = nf$nightPassage_rates, wc = wc,
										  wc_prop = 5/6, stratAssign_comp = stratComp,
										  stratAssign_night = stratNight, boots = numBoots)

if(spp == "STHD"){
	trap <- trap %>% mutate(lenCat = ifelse(LGDFLmm >= 780, "Lg", "Sm"))
	w <- c("GenStock", "lenCat")
	h <- c("releaseGroup")
	hnc <- c("releaseGroup", "GenSex") # testing HNC
} else {
	trap <- trap %>% mutate(GenSex = ifelse(GenSex == "U", NA, GenSex))
	w <- c("GenStock", "GenSex")
	h <- c("releaseGroup")
	hnc <- c("releaseGroup")
}

est_comp <- list()
tagRates <- tagRates %>% filter(!is.na(tagRate), tagRate > 0)
print("old")
set.seed(7)
est_comp[["old"]] <- ascension_composition(trap = trap, stratAssign_comp = stratComp, boots = numBoots,
														 pbt_var = "releaseGroup", tagRates = tagRates,
														 H_vars = h,
														 HNC_vars = hnc,
														 W_vars = w, wc_binom = exp_wc)
print("Acc")
set.seed(7)
est_comp[["Acc"]] <- HNC_expand(trap = trap, stratAssign_comp = stratComp, boots = numBoots,
										  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
										  adclip_var = "LGDMarkAD", tagRates = tagRates,
										  H_vars = h,
										  HNC_vars = hnc,
										  W_vars = w, wc_binom = exp_wc, method = "Acc")
print("MLE")
set.seed(7)
est_comp[["MLE"]] <- HNC_expand(trap = trap, stratAssign_comp = stratComp, boots = numBoots,
										  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
										  adclip_var = "LGDMarkAD", tagRates = tagRates,
										  H_vars = h,
										  HNC_vars = hnc,
										  W_vars = w, wc_binom = exp_wc, method = "MLE")

print("MLE GSI variable")
set.seed(7)
# randomizing GSI draws so point estimate draws are not correlated
gsiDraws <- gsiDraws %>% select(MasterID, sample(2:ncol(.), size = ncol(.) - 1, replace = FALSE))
set.seed(7)
est_comp[["MLE_gsiVar"]] <- HNC_expand_unkGSI(trap = trap, stratAssign_comp = stratComp, boots = numBoots,
															 pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
															 adclip_var = "LGDMarkAD",
															 sampID = "MasterID",
															 tagRates = tagRates,
															 H_vars = h,
															 HNC_vars = hnc,
															 W_vars = w, wc_binom = exp_wc,
															 GSI_draws = gsiDraws, n_point = ceiling(.1 * numBoots), GSI_var = "GenStock",
															 method = "MLE")

est_comp[["MLE"]][[1]] %>% filter(rear == "HNC") %>% full_join(
	est_comp[["Acc"]][[1]] %>% filter(rear == "HNC") %>% rename(acc = total)
) %>% filter(abs(total - acc) > 1)



est_escp <- list()
for(i in 1:4){
	if(exists("sthd_splitPITdata")) est_comp[[i]] <- splitByPIT(est_comp[["MLE_gsiVar"]], sthd_splitPITdata)

	templates <- apply_fallback_rates(breakdown = est_comp[[i]], fallback_rates = nf$fallback_rates,
												 split_H_fallback = "var1",
												 split_HNC_fallback = "var1",
												 split_W_fallback = "var1",
												 H_groups = NULL, HNC_groups = NULL, W_groups = NULL,
												 stratAssign_fallback = stratFallback, stratAssign_comp = stratComp, alpha_ci = .1,
												 output_type = "summary")
	if(spp == "STHD"){
		# assign upper and lower stock group
		templates$W$stockGroup <- ifelse(templates$W$var1 == "LSNAKE", "lower", "upper")
		templates$HNC <- templates$HNC %>% select(-stockGroup) %>% left_join(sthd_stockGroup, by = c("var1" = "releaseGroup"))
		templates$HNC <- templates$HNC %>% mutate(stockGroup = ifelse(is.na(stockGroup) & grepl("_lower$", var1), "lower",
																						  ifelse(is.na(stockGroup) & grepl("_upper$", var1), "upper", stockGroup)))
		templates$H <- templates$H %>% select(-stockGroup) %>% left_join(sthd_stockGroup, by = c("var1" = "releaseGroup"))
		templates$H <- templates$H %>% mutate(stockGroup = ifelse(is.na(stockGroup) & grepl("_lower$", var1), "lower",
																					 ifelse(is.na(stockGroup) & grepl("_upper$", var1), "upper", stockGroup)))
	} else {
		templates$W$stockGroup <- "upper"
		templates$HNC$stockGroup <- "upper"
		templates$H$stockGroup <- "upper"
	}

	# all should have 0 rows
	if((templates$HNC %>% filter(is.na(stockGroup)) %>% nrow) != 0) stop("HNC")
	if((templates$H %>% filter(is.na(stockGroup)) %>% nrow) != 0) stop("H")
	if((templates$W %>% filter(is.na(stockGroup)) %>% nrow) != 0) stop("W")

	est_escp[[names(est_comp)[i]]] <- apply_fallback_rates(breakdown = est_comp[[i]], fallback_rates = nf$fallback_rates,
																			 split_H_fallback = "var1",
																			 split_HNC_fallback = "var1",
																			 split_W_fallback = "var1",
																			 H_groups = templates$H, HNC_groups = templates$HNC, W_groups = templates$W,
																			 stratAssign_fallback = stratFallback, stratAssign_comp = stratComp, alpha_ci = .1,
																			 output_type = "summary")
}




load("../data/rda/STHD_sY19.rda")
usethis::use_data(fullReComplete, fullNi, stratFallback, stratNight, wc, stratComp,
						tagRates, trap, gsiDraws, sthd_splitPITdata, sthd_stockGroup)
