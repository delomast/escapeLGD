library(tidyverse)
library(lubridate)
load("example_data_STHD18.rda")

nf <- nightFall(wc = wc, full_reascend = full_reascend_flip, full_night = full_night,
					 stratAssign_fallback = stratAssign, stratAssign_night = stratAssign,
					 full_spillway = NULL, boots = 20)
exp_wc <- expand_wc_binom_night(nightPassage_rates = nf$nightPassage_rates, wc = wc,
										  wc_prop = 5/6, stratAssign_comp = stratAssign_comp,
										  stratAssign_night = stratAssign, boots = 20)




est_comp <- ascension_composition(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
											 pbt_var = "releaseGroup", tagRates = tagRates,
											 H_vars = c("releaseGroup", "GenSex"),
											 HNC_vars = c("releaseGroup"),
											 W_vars = c("GenStock", "GenSex"), wc_binom = exp_wc)

# test one strata function
HNC_expand_one_strat_MLE(trap = trap, pbt_var = "releaseGroup", tagRates = bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1)),
							H_vars = c("releaseGroup", "GenSex"),
							HNC_vars = c("releaseGroup"),
							W_vars = c("GenStock", "GenSex"), wc_expanded = 10000
)
HNC_expand_one_strat(trap = trap, pbt_var = "releaseGroup", tagRates = bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1)),
							H_vars = c("releaseGroup", "GenSex"),
							HNC_vars = c("releaseGroup"),
							W_vars = c("GenStock", "GenSex"), wc_expanded = 10000
)

debugonce(HNC_expand_one_strat)
HNC_expand_one_strat_MLE(trap = trap %>% filter(sWeek %in% (stratAssign_comp %>% filter(stratum == 1) %>% pull(sWeek)),
														  releaseGroup != "Unassigned"),
							pbt_var = "releaseGroup", tagRates = bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1)),
							H_vars = c("releaseGroup", "GenSex"),
							HNC_vars = c("releaseGroup"),
							W_vars = c("GenStock", "GenSex"), wc_expanded = 10000
)$estimates %>% filter(is.na(var2)) %>% as.data.frame() %>% mutate(total = round(total, 3)) %>%
	group_by(var1) %>% summarise(total = sum(total)) %>%
	full_join(trap %>% filter(sWeek %in% (stratAssign_comp %>% filter(stratum == 1) %>% pull(sWeek)),
									  releaseGroup != "Unassigned") %>% count(releaseGroup) %>% rename(var1 = releaseGroup)) %>% as.data.frame()



HNC_expand_one_strat(trap = trap %>% filter(sWeek %in% (stratAssign_comp %>% filter(stratum == 1) %>% pull(sWeek)),
														  releaseGroup != "Unassigned"),
							pbt_var = "releaseGroup", tagRates = bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1)),
							H_vars = c("releaseGroup", "GenSex"),
							HNC_vars = c("releaseGroup"),
							W_vars = c("GenStock", "GenSex"), wc_expanded = 10000
)$estimates %>% filter(is.na(total))


est_comp <- HNC_expand(trap = trap, stratAssign_comp = stratAssign_comp, boots = 20,
							  pbt_var = "releaseGroup", timestep_var = "sWeek", physTag_var = "physTag",
							  adclip_var = "LGDMarkAD", tagRates = tagRates,
			  H_vars = c("releaseGroup", "GenSex"),
			  HNC_vars = c("releaseGroup"),
			  W_vars = c("GenStock", "GenSex"), wc_binom = exp_wc)
sum(is.na(est_comp[[1]]$total))
est_comp[[1]] %>% filter(is.na(total))

sum(is.na(est_comp[[2]]$total))
est_comp[[2]] %>% filter(is.na(total)) %>% pull(total) %>% unique
est_comp[[1]] %>% filter(is.na(total)) %>% pull(total) %>% unique

est_comp[[2]] %>% filter(is.na(total)) %>% count(rear, var1)
est_comp[[1]] %>% filter(is.na(total)) %>% count(rear, var1)

trap %>% filter(sWeek %in% (stratAssign_comp %>% filter(stratum == 1) %>% pull(sWeek))) %>%
	filter(LGDMarkAD == "AD") %>% count(releaseGroup) %>% as.data.frame()


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

