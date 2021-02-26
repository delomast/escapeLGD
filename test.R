library(tidyverse)
library(lubridate)
load("example_data_STHD18.rda")

nf <- nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)
nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)

# fake_spillway <- full_reascend_flip %>% mutate(totalFall = 5, laterAscend = 4) %>% select(sWeek, stockGroup, totalFall, laterAscend)
# nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
# 			 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
# 			 full_spillway = fake_spillway, boots = 20)


exp_wc <- expand_wc_binom_night(nightPassage_rates = nf$nightPassage_rates, wc = wc,
										  wc_prop = 5/6, stratAssign_comp = stratAssign_comp,
										  stratAssign_night = stratAssign, boots = 20)




est_comp <- ascension_composition(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
											 pbt_var = "releaseGroup", tagRates = tagRates,
											 H_vars = c("releaseGroup", "GenSex"),
											 HNC_vars = c("releaseGroup"),
											 W_vars = c("GenStock", "GenSex"), wc_binom = exp_wc)

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
			  W_vars = c("GenSex", "lenCat"), wc_binom = exp_wc, method = "MLE", testing = "manual")

est_comp[[1]] %>% rename(acTotal = total) %>% full_join(est_comp_MLE[[1]]) %>%
	mutate(total = round(total, 4), diff = acTotal - total) %>%
	filter(abs(diff) > 1) %>% arrange(var1) %>% as.data.frame() %>%
	filter(rear == "W") %>% arrange(stratum)

est_comp[[1]] %>% rename(acTotal = total) %>% full_join(est_comp_MLE[[1]]) %>%
	mutate(total = round(total, 4), diff = acTotal - total) %>%
	filter(rear == "W") %>% arrange(stratum, var1, var2) %>% as.data.frame() %>%
	filter(!is.na(var2))


est_comp[[1]] %>% filter(rear == "HNC", is.na(var2)) %>% pull(total) %>% sum
est_comp_MLE[[1]] %>% filter(rear == "HNC", is.na(var2)) %>% pull(total) %>% sum

est_comp[[1]] %>% filter(rear == "H", !is.na(var2)) %>% pull(total) %>% sum
est_comp_MLE[[1]] %>% filter(rear == "H", !is.na(var2)) %>% pull(total) %>% sum

est_comp[[1]] %>% filter(is.na(var2)) %>% pull(total) %>% sum
est_comp_MLE[[1]] %>% filter(is.na(var2)) %>% pull(total) %>% sum


est_comp[[1]] %>% filter(rear == "W", !is.na(var2)) %>% pull(total) %>% sum
est_comp_MLE[[1]] %>% filter(rear == "W", !is.na(var2)) %>% pull(total) %>% sum



trap %>% left_join(stratAssign_comp) %>% filter(stratum == 1) %>% filter(!is.na(pbtAssign), !pbtAssign, !physTag, LGDMarkAD == "AI") %>%
	count(GenSex, lenCat)




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
est_escp
est_comp[[1]] %>% filter(rear == "H", !is.na(var2)) %>% as.data.frame()

