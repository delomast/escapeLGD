library(tidyverse)
library(lubridate)
load("example_data_STHD18.rda")

nf <- nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)
nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)

fake_spillway <- full_reascend_flip %>% mutate(totalFall = 5, laterAscend = 4) %>% select(sWeek, stockGroup, totalFall, laterAscend)
nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
			 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
			 full_spillway = fake_spillway, boots = 20)


exp_wc <- expand_wc_binom_night(nightPassage_rates = nf$nightPassage_rates, wc = wc,
										  wc_prop = 5/6, stratAssign_comp = stratAssign_comp,
										  stratAssign_night = stratAssign, boots = 20)




est_comp <- ascension_composition(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
											 pbt_var = "releaseGroup", tagRates = tagRates,
											 H_vars = c("releaseGroup", "GenSex"),
											 HNC_vars = c("releaseGroup"),
											 W_vars = c("GenStock", "GenSex"), wc_binom = exp_wc)

# makeing categorical variable for testing
trap <- trap %>% mutate(lenCat = ifelse(LGDFLmm > 650, "large", ifelse(LGDFLmm > 600, "medium", "small")))

system.time(
est_comp <- HNC_expand(trap = trap, stratAssign_comp = stratAssign_comp, boots = 0,
							  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
							  adclip_var = "LGDMarkAD", tagRates = tagRates,
			  H_vars = c("releaseGroup", "GenSex"),
			  HNC_vars = c("releaseGroup"),
			  W_vars = c("GenSex", "lenCat"), wc_binom = exp_wc, method = "MLE", testing = "manual")
)

templates <- apply_fallback_rates(breakdown = est_comp, fallback_rates = nf$fallback_rates,
											 split_H_fallback = "var2",
											 split_HNC_fallback = "var1",
											 split_W_fallback = "both",
											 H_groups = NULL, HNC_groups = NULL, W_groups = NULL,
											 stratAssign_fallback, stratAssign_comp, alpha_ci = .1,
											 output_type = "summary")

for(i in 1:3) templates[[i]]$stockGroup <- sample(unique(nf$fallback_rates[[1]]$stockGroup), nrow(templates[[i]]), replace = TRUE)
# boots = 20
est_escp <- apply_fallback_rates(breakdown = est_comp, fallback_rates = nf$fallback_rates,
											split_H_fallback = "var2",
											split_HNC_fallback = "var1",
											split_W_fallback = "both",
											H_groups = templates$H, HNC_groups = templates$HNC, W_groups = templates$W,
											stratAssign_fallback = stratAssign, stratAssign_comp = stratAssign_comp, alpha_ci = .1,
											output_type = "summary")

options(error = recover)
options(error = NULL)



