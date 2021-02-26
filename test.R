library(tidyverse)
library(lubridate)
load("example_data_STHD18.rda")

nf <- nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)

exp_wc <- expand_wc_binom_night(nightPassage_rates = nf$nightPassage_rates, wc = wc,
										  wc_prop = 5/6, stratAssign_comp = stratAssign_comp,
										  stratAssign_night = stratAssign, boots = 20)

est_comp_old <- ascension_composition(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
											 pbt_var = "releaseGroup", tagRates = tagRates,
											 H_vars = c("releaseGroup", "GenSex"),
											 HNC_vars = c("releaseGroup"),
											 W_vars = c("GenSex", "lenCat"), wc_binom = exp_wc)

# making categorical variable for testing
trap <- trap %>% mutate(lenCat = ifelse(LGDFLmm > 650, "large", ifelse(LGDFLmm > 600, "medium", "small")))

est_comp <- HNC_expand(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
							  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
							  adclip_var = "LGDMarkAD", tagRates = tagRates,
			  H_vars = c("releaseGroup", "GenSex"),
			  HNC_vars = c("releaseGroup"),
			  W_vars = c("GenSex", "lenCat"), wc_binom = exp_wc, method = "Acc", testing = "manual")


est_comp_MLE <- HNC_expand(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
							  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
							  adclip_var = "LGDMarkAD", tagRates = tagRates,
			  H_vars = c("releaseGroup", "GenSex"),
			  HNC_vars = c("releaseGroup"),
			  W_vars = c("GenSex", "lenCat"), wc_binom = exp_wc, method = "MLE", testing = "complete")


templates <- apply_fallback_rates(breakdown = est_comp_MLE, fallback_rates = nf$fallback_rates,
											 split_H_fallback = "var2",
											 split_HNC_fallback = "var1",
											 split_W_fallback = "both",
											 H_groups = NULL, HNC_groups = NULL, W_groups = NULL,
											 stratAssign_fallback, stratAssign_comp, alpha_ci = .1,
											 output_type = "summary")

for(i in 1:3) templates[[i]]$stockGroup <- sample(unique(nf$fallback_rates[[1]]$stockGroup), nrow(templates[[i]]), replace = TRUE)
# boots = 20
est_escp <- apply_fallback_rates(breakdown = est_comp_MLE, fallback_rates = nf$fallback_rates,
											split_H_fallback = "var2",
											split_HNC_fallback = "var1",
											split_W_fallback = "both",
											H_groups = templates$H, HNC_groups = templates$HNC, W_groups = templates$W,
											stratAssign_fallback = stratAssign, stratAssign_comp = stratAssign_comp, alpha_ci = .1,
											output_type = "summary")

options(error = recover)
options(error = NULL)
est_escp


trap <- trap %>% mutate(fakeGSI = sample(c("g1", "g2"), nrow(.), prob = c(1, 2), replace = TRUE),
								Ind = paste0("Ind_", 1:nrow(.)))

gsiDraws <- trap %>% select(Ind)
for(i in 1:2000){
	gsiDraws <- gsiDraws %>% bind_cols(sample(c("g1", "g2"), nrow(.), prob = c(1, 2), replace = TRUE)
	)
}

gsiDraws[[3]] <- "g3"

est_comp_t1 <- HNC_expand_unkGSI(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
									pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
									adclip_var = "LGDMarkAD", sampID = "Ind",
									tagRates = tagRates,
									H_vars = c("releaseGroup", "GenSex"),
									HNC_vars = c("releaseGroup"),
									W_vars = c("GenSex", "fakeGSI"), wc_binom = exp_wc,
									GSI_draws = gsiDraws, n_point = 1, GSI_var = "fakeGSI",
									method = "MLE")

est_comp_t1[[1]] %>% filter(!is.na(var2), rear == "W")
est_comp_t1[[1]] %>% filter(is.na(var2), rear == "W")




templates <- apply_fallback_rates(breakdown = est_comp_t1, fallback_rates = nf$fallback_rates,
											 split_H_fallback = "var2",
											 split_HNC_fallback = "var1",
											 split_W_fallback = "both",
											 H_groups = NULL, HNC_groups = NULL, W_groups = NULL,
											 stratAssign_fallback, stratAssign_comp, alpha_ci = .1,
											 output_type = "summary")
for(i in 1:3) templates[[i]]$stockGroup <- sample(unique(nf$fallback_rates[[1]]$stockGroup), nrow(templates[[i]]), replace = TRUE)
est_escp_uG <- apply_fallback_rates(breakdown = est_comp_t1, fallback_rates = nf$fallback_rates,
											split_H_fallback = "var2",
											split_HNC_fallback = "var1",
											split_W_fallback = "both",
											H_groups = templates$H, HNC_groups = templates$HNC, W_groups = templates$W,
											stratAssign_fallback = stratAssign, stratAssign_comp = stratAssign_comp, alpha_ci = .1,
											output_type = "summary")
