# functions used in estimating composition with the trap biological data

#' estimates composition for a non-PBT categorical variable
#' @param values values of the category for all the samples, a vector
#' @param boots number of bootstrap iterations
#' @return a list with two components: first is tibble of point estimates,
#'   second is matrix with bootstraps. Columns are categories are in order following the
#'   point estimates and are named accordingly. Rows are bootstrap iterations.
#' @keywords internal
#' @noRd
nonPBT_breakdown <- function(values, boots){
	props <- table(values) # removes NAs automatically
	nam <- names(props)
	sampleSize <- sum(props)
	props <- as.numeric(props/sampleSize)
	output <- list(tibble(group = nam, prop = props))
	output[[2]] <- t(rmultinom(n = boots, size = sampleSize, prob = props)) / sampleSize
	colnames(output[[2]]) <- nam
	return(output)
}

# correcting W for HNC untagged
# this should be optional if easy to do
# to break down HNC, can pass phystag variable and total untagged count - should
# have bootstraps from H,HNC,W breakdown. Then expand appropriately splitting
# phytag and no phystag

# then for wild pass the HNC estimates and subtract as appropriate
# wild composition just has to be estimated after HNC


# also incorporate GSI uncertainty



#' estimates composition for PBT groups
#' @param values values of the category for all the samples, a vector
#' @param tagRates a tibble with release group in col1 and tag rate in col2
#' @return returns tibble with group (release Group or Unassigned) and
#'   prop (proportion of the stratum) columns
#' @keywords internal
#' @noRd
PBT_expand_calc <- function(values, tagRates){
	props <- table(values) # removes NAs automatically
	pr <- tibble(group = names(props), count = as.numeric(props)) %>% left_join(tagRates, by = "group") %>%
		mutate(expand = count / tagRate, diff = expand - count)
	if(!"Unassigned" %in% pr$group) pr <- pr %>% bind_rows(tibble(group = "Unassigned",
																					  count = 0, tagRate = 1, expand = 0, diff = 0))
	# truncate expansions when not enough Unassigned
	if((pr %>% filter(group == "Unassigned") %>% pull(count)) < sum(pr$diff)){
		pr$diff <- (pr %>% filter(group == "Unassigned") %>% pull(count)) * (pr$diff / sum(pr$diff))
		pr$expand <- pr$count + pr$diff
		pr$expand[pr$group == "Unassigned"] <- 0
	} else {
		pr$expand[pr$group == "Unassigned"] <- pr$expand[pr$group == "Unassigned"] - sum(pr$diff)
	}
	return(pr %>% mutate(prop = expand / sum(expand)) %>% select(group, prop))
}

#' estimates composition for a non-PBT categorical variable
#' @param values values of the category for all the samples, a vector
#' @param tagRates a tibble with release group in col1 and tag rate in col2
#' @param boots number of bootstrap iterations
#' @return a list with two components: first is tibble of point estimates,
#'   second is matrix with bootstraps. Columns are categories are in order following the
#'   point estimates and are named accordingly. Rows are bootstrap iterations.
#' @keywords internal
#' @noRd
PBT_breakdown <- function(values, tagRates, boots){
	colnames(tagRates) <- c("group", "tagRate")
	if("Unassigned" %in% tagRates[[1]]){
		if(tagRates[[2]][tagRates[[1]] == "Unassigned"] != 1){
			stop("Unassigned must not be included in the tag rate file or have a tag rate of 1")
		}
	} else {
		# add Unassigned with tag rate of 1
		tagRates <- bind_rows(tagRates, tibble(group = "Unassigned", tagRate = 1))
	}
	output <- list(PBT_expand_calc(values, tagRates))
	# non-parametric bootstrapping
	output[[2]] <- matrix(nrow = boots, ncol = nrow(output[[1]]))
	for(i in 1:boots){
		temp <- PBT_expand_calc(values = sample(values, size = length(values), replace = TRUE), tagRates)
		output[[2]][i,] <- output[[1]] %>% select(group) %>% left_join(temp, by = "group") %>%
			mutate(prop = replace_na(prop, 0)) %>% pull(prop)
	}
	colnames(output[[2]]) <- output[[1]]$group
	return(output)
}

#' estimates composition of one stratum for one or two variables
#' @param trapStratumData values of the category for all the samples, a vector
#' @param vars a tibble with release group in col1 and tag rate in col2
#' @param pbt_var name of pbt variable (character)
#' @return a list with one component for each variable.They contain point
#'   estimates and bootstrap estimates
#' @keywords internal
#' @noRd
subGroup_breakdown <- function(trapStratumData, vars, pbt_var, tagRates, boots = boots){
	list_break <- list()
	v1 <- vars[1]
	v2 <- if(length(vars) == 2) vars[2] else NULL
	trapStratumData <- trapStratumData[!is.na(trapStratumData[[v1]]),]
	if(v1 == pbt_var){
		list_break[[1]] <- PBT_breakdown(values = trapStratumData[[v1]], tagRates = tagRates, boots = boots)
	} else {
		list_break[[1]] <- nonPBT_breakdown(values = trapStratumData[[v1]], boots = boots)
	}
	names(list_break)[1] <- v1
	if(!is.null(v2)){
		trapStratumData <- trapStratumData[!is.na(trapStratumData[[v2]]),] # remove observations missing data for v2
		# storing all the values for the subvariables
		# as a massive list and named/indexed according to order
		list_break[[2]] <- list() # index here is category of the first variable
		for(j in 1:nrow(list_break[[1]][[1]])){ # run breakdown for all categories
			v1Cat <- list_break[[1]][[1]]$group[j] # current value of v1
			tempData <- trapStratumData[[v2]][trapStratumData[[v1]] == v1Cat]
			if(v2 == pbt_var){
				list_break[[2]][[j]] <- PBT_breakdown(values = tempData, tagRates = tagRates, boots = boots)
			} else {
				list_break[[2]][[j]] <- nonPBT_breakdown(values = tempData, boots = boots)
			}
		}
		names(list_break[[2]]) <- list_break[[1]][[1]]$group
		names(list_break)[2] <- v2
	}
	return(list_break)
}
