# functions for fallback adn reascension estimation

#' log-likelihood of fallback and reascension rates
#' @param pf P(fallback)
#' @param pre_f P(reascend|fallback)
#' @param dfr count that reascended after known fallback (spillway)
#' @param df count that fellback (spillway)
#' @param dr count of reascensions (ladder)
#' @param dt count of ascensions (ladder)
#' @keywords internal
#' @noRd
#' @return log likelihood value
fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt){
	return(dbinom(dfr, df, pre_f, log = TRUE) + dbinom(dr, dt, pf*pre_f, log = TRUE))
}

#' gradient function of fallback_log_likelihood
#' @param par c(P(fallback), P(reascend|fallback))
#' @param dfr count that reascended after known fallback (spillway)
#' @param df count that fellback (spillway)
#' @param dr count of reascensions (ladder)
#' @param dt count of ascensions (ladder)
#' @keywords internal
#' @noRd
#' @return gradient
gradient_fallback_log_likelihood <- function(par, dfr, df, dr, dt){
	pf <- par[1]
	pre_f <- par[2]

	# try(
	# print(data.frame(an = c(
	# 	# partial deriv with respect to pf
	# 	(dr / pf) - (((dt - dr) * pre_f) / (1 - (pf * pre_f))),
	# 	# partial deriv with respect to pre_f
	# 	(dfr / pre_f) - ((df - dfr)/(1 - pre_f)) + (dr / pre_f) -
	# 		(((dt - dr) * pf) / (1 - (pf*pre_f)))
	# ),
	# 					  num = numDeriv::grad(optimllh, par, dfr = dfr, df = df, dr = dr, dt = dt)) %>%
	# 					  	mutate(diff = an - num)
	# 					  )
	# )

	return(c(
		# partial deriv with respect to pf
		(dr / pf) - (((dt - dr) * pre_f) / (1 - (pf * pre_f))),
		# partial deriv with respect to pre_f
		(dfr / pre_f) - ((df - dfr)/(1 - pre_f)) + (dr / pre_f) -
			(((dt - dr) * pf) / (1 - (pf*pre_f)))
	))
}

#' wrapper for fallback_log_likelihood that keeps all variables in par
#' @param par c(P(fallback), P(reascend|fallback))
#' @param dfr count that reascended after known fallback (spillway)
#' @param df count that fellback (spillway)
#' @param dr count of reascensions (ladder)
#' @param dt count of ascensions (ladder)
#' @keywords internal
#' @noRd
#' @return log likelihood value
optimllh <- function(par, dfr, df, dr, dt){
	# initial testing with softmax and NM or BFGS w/o gradient
	# pf <- exp(par[1])
	# pre_f <- exp(par[2])
	# pf <- pf / (pf + 2.718282)
	# pre_f <- pre_f / (pre_f + 2.718282)
	pf <- par[1]
	pre_f <- par[2]
	return(fallback_log_likelihood(pf, pre_f, dfr, df, dr, dt))
}

