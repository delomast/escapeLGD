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


#' softmax function with numerical stability
#' Values will sum to 1
#' may return values of 0, but not Inf
#' @param x real numbers to transform into simplex via softmax
#' @keywords internal
#' @noRd
softMax <- function(x){
	pr <- exp(x - max(x))
	return(pr / sum(pr))
}


#' MLE estimates composition for PBT groups
#' @param values values of the category for all the samples, a vector
#' @param tagRates a tibble with release group as "group" and tag rate as "tagRate"
#' @return returns tibble with group (release Group or Unassigned) and
#'   prop (proportion of the stratum) columns
#' @keywords internal
#' @noRd
PBT_expand_calc_MLE <- function(values, tagRates){
	values <- values[!is.na(values)]
	numUnassign <- sum(values == "Unassigned")
	numAssign <- table(values[values != "Unassigned"])
	# no PBT assigned fish
	if(length(numAssign) < 1) return(tibble(group = "Unassigned", prop = 1))

	pr <- tibble(group = names(numAssign), count = as.numeric(numAssign)) %>% left_join(tagRates, by = "group")

	# now MLE estimation of proportions\
	propsMLE <- optim(par = rep(1, nrow(pr) + 1), fn = PBT_optimllh, gr = PBT_grad,
								nGroups = pr$count, nUntag = numUnassign, tagRates = pr$tagRate,
							  control = list(fnscale = -1), method = "BFGS"
							)
	if(propsMLE$convergence != 0) warning("Convergence error in PBT_expand_calc_MLE")
	# apply softmax to translate into probabilities
	par_prob <- softMax(propsMLE$par)
	# fixing computational limit on inferring parameters to be zero
	if(numUnassign == 0){
		par_prob[length(par_prob)] <- 0
		par_prob <- par_prob / sum(par_prob)
	}
	return(tibble(group = c(pr$group, "Unassigned"), prop = par_prob))
}

#' log-likelihood of PBT groups
#' @param pGroups proportion belonging to each PBT group
#' @param pW proportion that is wild
#' @param nGroups counts for each PBT group
#' @param nUntag counts of untagged
#' @param tagRates tagRates in order of pGroups
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_log_likelihood <- function(pGroups, pW, nGroups, nUntag, tagRates){
	y <- pW + sum((1 - tagRates) * pGroups)
	pGroups_tag <- pGroups * tagRates
	return(dmultinom(x = c(nGroups, nUntag), prob = c(pGroups_tag, y), log = TRUE))
}

#' wrapper for PBT log-likelihood that keeps all variables in par
#' This requires at least one PBT tagged group is present
#' @param par real numbers corresponding to softmax probabilities of c(PBT groups, wild)
#' @param nGroups counts for each PBT group
#' @param nUntag counts of untagged
#' @param tagRates tagRates in order of pGroups
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_optimllh <- function(par, nGroups, nUntag, tagRates){
	# apply softmax to translate into probabilities
	par_prob <- softMax(par)
	pGroups <- par_prob[1:(length(par_prob) - 1)]
	pW <- par_prob[length(par_prob)]
	return(PBT_log_likelihood(pGroups = pGroups, pW = pW,
									  nGroups = nGroups, nUntag = nUntag, tagRates = tagRates))
}

PBT_grad <- function(par, nGroups, nUntag, tagRates){
	# apply softmax to translate into probabilities
	par_prob <- softMax(par)
	nPBT <- length(par_prob) - 1
	pGroups <- par_prob[1:nPBT]
	pW <- par_prob[nPBT + 1]
	# calculate probs for multinomial
	y <- pW + sum((1 - tagRates) * pGroups)

	gr <- rep(0, length(par))

	for(i in 1:(nPBT + 1)){
		smDeriv <- -par_prob * par_prob[i]
		smDeriv[i] <- par_prob[i] * (1 - par_prob[i])
		gr[i] <- sum((nGroups / pGroups) * smDeriv[1:nPBT]) +
			((nUntag / y) * (sum((1 - tagRates[1:nPBT]) * smDeriv[1:nPBT]) + smDeriv[nPBT + 1]))
	}

	# print(data.frame(an = gr,
	# 					  num = numDeriv::grad(PBT_optimllh, par, nGroups = nGroups, nUntag = nUntag, tagRates = tagRates)) %>%
	# 					  	mutate(diff = an - num)
	# 					  )

	return(gr)
}



#' MLE estimates composition for var2 within PBT groups
#' @param tagRates a tibble with release group as "group" and tag rate as "tagRate"
#' @return returns a matrix with PBT group/Unassigned as rows and var2 as cols, and props within each group as values
#' @keywords internal
#' @noRd
PBT_var2_calc_MLE <- function(v2Data, piGroup, tagRates){
	# now MLE estimation of proportions\
	propsMLE <- optim(par = rep(1, nrow(v2Data) * ncol(v2Data)), fn = PBT_var2_optimllh, gr = PBT_var2_grad,
							piGroup = piGroup, v2Data = v2Data, tagRates = tagRates,
							control = list(fnscale = -1), method = "BFGS"
					)
	if(propsMLE$convergence != 0) warning("Convergence error in PBT_var2_calc_MLE")
	# apply softmax to translate into probabilities
	varProbs <- matrix(propsMLE$par, nrow = nrow(v2Data), ncol = ncol(v2Data), byrow = TRUE)
	for(i in 1:nrow(varProbs)){
		varProbs[i,] <- softMax(varProbs[i,])
	}
	rownames(varProbs) <- rownames(v2Data)
	colnames(varProbs) <- colnames(v2Data)
	return(varProbs)
}

#' wrapper for estimating var 2 with PBT as var1
#' @param par real numbers corresponding to softmax probabilities of
#' @param piGroup "known" proportions of each group in order of rows of v2Data
#' @param v2Data counts of each combination
#' @param tagRates tagRates in order of piGroup
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_var2_optimllh <- function(par, piGroup, v2Data, tagRates){
	# apply softmax to translate into probabilities
	varProbs <- matrix(par, nrow = nrow(v2Data), ncol = ncol(v2Data), byrow = TRUE)
	for(i in 1:nrow(varProbs)){
		varProbs[i,] <- softMax(varProbs[i,])
	}
	return(PBT_var2_log_likelihood(varProbs = varProbs, piGroup = piGroup,
											 v2Data = v2Data, tagRates = tagRates))
}

#' log_likelihood for estimating var 2 with PBT as var1
#' @param varProbs rows are groups, cols are var2 cats, values are probs within each group
#' @param piGroup "known" proportions of each group in order of rows of v2Data
#' @param v2Data counts of each combination
#' @param tagRates tagRates in order of piGroup
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_var2_log_likelihood <- function(varProbs, piGroup,
												v2Data, tagRates){
	# product of multinomial likelihoods (sum of log-likelihoods)
	llh <- 0
	nPBT <- nrow(v2Data) - 1
	# first the tagged fish
	for(i in 1:nPBT){
		llh <- llh + dmultinom(x = v2Data[i,], prob = varProbs[i,], log = TRUE)
	}
	# now the untagged fish
	untagProb <- piGroup
	untagProb[1:nPBT] <- untagProb[1:nPBT] * (1 - tagRates[1:nPBT])
	untagProb <- untagProb / sum(untagProb)
	# untagVarProbs <- matrix(nrow = nrow(varProbs), ncol = ncol(varProbs))
	# for(i in 1:nrow(v2Data)) untagVarProbs[i,] <- untagProb[i] * varProbs[i,]
	# untagVarProbs <- colSums(untagVarProbs)
	untagVarProbs <- drop(matrix(untagProb, nrow = 1) %*% varProbs)
	llh <- llh + dmultinom(x = v2Data[nrow(v2Data),], prob = untagVarProbs, log = TRUE)

	return(llh)
}

PBT_var2_grad <- function(par, piGroup, v2Data, tagRates){

	gr <- rep(0, length(par))
	nPBT <- nrow(v2Data) - 1

	# apply softmax to translate into probabilities
	varProbs <- matrix(par, nrow = nrow(v2Data), ncol = ncol(v2Data), byrow = TRUE)
	for(i in 1:nrow(varProbs)){
		varProbs[i,] <- softMax(varProbs[i,])
	}

	# first the tagged fish
	pos <- 1
	for(i in 1:nPBT){
		for(j in 1:ncol(v2Data)){
			gr[pos] <- sum(v2Data[i,] * (((1:ncol(v2Data)) == j) - varProbs[i,j]))
			pos <- pos + 1
		}
	}
	# gr[1:(nPBT * ncol(v2Data))] <- as.vector(t(v2Data[1:nPBT,] * (1 - varProbs[1:nPBT,])))

	# now the untagged fish
	untagProb <- piGroup
	untagProb[1:nPBT] <- untagProb[1:nPBT] * (1 - tagRates[1:nPBT])
	untagProb <- untagProb / sum(untagProb)

	untagVarProbs <- drop(matrix(untagProb, nrow = 1) %*% varProbs)

	unTagData <- v2Data[nPBT + 1,]
	pos <- 1
	for(p in 1:(nPBT+1)){
		for(j in 1:ncol(v2Data)){
			smDeriv <- varProbs[p,] * (((1:ncol(v2Data)) == j) - varProbs[p,j])
			gr[pos] <- gr[pos] + sum((unTagData / untagVarProbs) * untagProb[p] * smDeriv)
			pos <- pos + 1
		}
	}

	# print(data.frame(an = gr,
	# 					  num = numDeriv::grad(PBT_var2_optimllh, par, piGroup = piGroup, v2Data = v2Data, tagRates = tagRates)) %>%
	# 					  	mutate(diff = an - num) %>% pull(diff) %>% mean
	# 					  )

	return(gr)
}







#' MLE estimates composition for var3 within W
#' @param piGroup "known" proportions of each group in order of rows of
#' @param tagRates tagRates in order of piGroup
#' @param var2Probs known composition of var2(of 3). rows are PBT groups with Unassigned (wild) at the end
#' @param hnc_var3 counts of PBT assigned in cats of var3
#' @param un_counts counts of unassigned in cats of var2 (rows) and var3 (cols)
#' @return returns a matrix with var2 as rows and var3 as cols, and props within each group as values
#' @keywords internal
#' @noRd
PBT_var3_calc_MLE <- function(piGroup, tagRates, var2Probs, hnc_var3, un_counts){
	# now MLE estimation of proportions\
	propsMLE <- optim(par = rep(1, (nrow(hnc_var3) * ncol(hnc_var3)) + (nrow(un_counts) * ncol(un_counts))),
							fn = PBT_var3_optimllh, # gr = , # need to add gradient
							piGroup = piGroup, var2Probs = var2Probs, tagRates = tagRates,
							hnc_var3 = hnc_var3, un_counts = un_counts,
							control = list(fnscale = -1), method = "BFGS"
	)
	if(propsMLE$convergence != 0) warning("Convergence error in PBT_var3_calc_MLE")
	# apply softmax to translate into probabilities
	var3Probs <- matrix(propsMLE$par, nrow = nrow(hnc_var3) + nrow(un_counts), ncol = ncol(hnc_var3), byrow = TRUE)
	# don't need HNC values
	var3Probs <- var3Probs[(nrow(hnc_var3) + 1):nrow(var3Probs),]
	for(i in 1:nrow(var3Probs)){
		var3Probs[i,] <- softMax(var3Probs[i,])
	}
	rownames(var3Probs) <- rownames(un_counts)
	colnames(var3Probs) <- colnames(un_counts)
	return(var3Probs)
}


#' wrapper for estimating var 3 in W
#' @param par real numbers corresponding to softmax probabilities of
#' @param piGroup "known" proportions of each group in order of rows of
#' @param tagRates tagRates in order of piGroup
#' @param var2Probs known composition of var2(of 3). rows are PBT groups with Unassigned (wild) at the end
#' @param hnc_var3 counts of PBT assigned in cats of var3
#' @param un_counts counts of unassigned in cats of var2 (rows) and var3 (cols)
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_var3_optimllh <- function(par, piGroup, tagRates, var2Probs, hnc_var3, un_counts){
	# apply softmax to translate into probabilities
	# these are the probabilities of var3 within the pbt groups and then
	# within the var2 groups of the unassigned fish
	# so, need one row for each PBT group and one for each var2 group of the Unassigned
	# and one column for each category of var3
	# PBT groups are first, and then cats of var2 in Unassigned
	var3Probs <- matrix(par, nrow = nrow(hnc_var3) + nrow(un_counts), ncol = ncol(hnc_var3), byrow = TRUE)
	for(i in 1:nrow(var3Probs)){
		var3Probs[i,] <- softMax(var3Probs[i,])
	}
	return(PBT_var3_log_likelihood(var3Probs = var3Probs, piGroup = piGroup,
											 tagRates = tagRates, var2Probs = var2Probs,
											 hnc_var3 = hnc_var3, un_counts = un_counts))
}


#' log_likelihood for estimating var 3 in W
#' @param var3Probs Probabilities of var3 cats (cols) in PBT groups and in var2 cats of W (rows)
#' @param piGroup "known" proportions of each group in order of rows of
#' @param tagRates tagRates in order of piGroup
#' @param var2Probs known composition of var2(of 3). rows are PBT groups with Unassigned (wild) at the end
#' @param hnc_var3 counts of PBT assigned in cats of var3
#' @param un_counts counts of unassigned in cats of var2 (rows) and var3 (cols)
#' @keywords internal
#' @noRd
#' @return log likelihood value
PBT_var3_log_likelihood <- function(var3Probs, piGroup, tagRates, var2Probs,
												hnc_var3, un_counts){
	# product of multinomial likelihoods (sum of log-likelihoods)
	llh <- 0
	nPBT <- nrow(hnc_var3)
	# first the tagged fish
	for(i in 1:nPBT){
		llh <- llh + dmultinom(x = hnc_var3[i,], prob = var3Probs[i,], log = TRUE)
	}
	# now the untagged fish

	untagProb <- piGroup
	untagProb[1:nPBT] <- untagProb[1:nPBT] * (1 - tagRates[1:nPBT])
	untagProb <- untagProb / sum(untagProb) # this is compostion of the Unassigned fish

	# this is probabilities an unassigned fish is in cat of var2 x var3
	untagVar3Probs <- matrix(0, nrow = nrow(un_counts), ncol = ncol(un_counts))
	for(i in 1:nrow(un_counts)){
		for(j in 1:nPBT){
			untagVar3Probs[i,] <- untagVar3Probs[i,] + (untagProb[j] * var2Probs[j,i] * var3Probs[j,])
		}
		untagVar3Probs[i,] <- untagVar3Probs[i,] + (untagProb[(nPBT+1)] * var2Probs[(nPBT+1),i] * var3Probs[(nPBT+i),])
	}

	llh <- llh + dmultinom(x = as.vector(un_counts), prob = as.vector(untagVar3Probs), log = TRUE)

	return(llh)
}





























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
