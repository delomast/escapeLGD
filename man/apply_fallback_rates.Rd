% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_fallback_rates_summarize.R
\name{apply_fallback_rates}
\alias{apply_fallback_rates}
\title{applies the estimated fallback rates to the estiamted ascension composition
to yield estimated escapement}
\usage{
apply_fallback_rates(
  breakdown,
  fallback_rates,
  split_H_fallback = c("var1", "var2", "both"),
  split_HNC_fallback = c("var1", "var2", "both"),
  split_W_fallback = c("var1", "var2", "both"),
  H_groups = NULL,
  HNC_groups = NULL,
  W_groups = NULL,
  stratAssign_fallback,
  stratAssign_comp,
  alpha_ci = 0.1,
  output_type = c("summary", "W_boot", "full")
)
}
\arguments{
\item{breakdown}{The output of either \code{HNC_expand_unkGSI}, \code{HNC_expand}, or
\code{ascension_composition}}

\item{fallback_rates}{The fallback rates portion of the output of \code{nightFall}}

\item{split_H_fallback}{one of "var1", "var2", or "both" indicating which variable(s)
should be used to link estimates of composition with fallback rates for "H" fish}

\item{H_groups}{\code{NULL} to get templates returned, or a tibble linking categories to stockGroups}

\item{HNC_groups}{\code{NULL} to get templates returned, or a tibble linking categories to stockGroups}

\item{W_groups}{\code{NULL} to get templates returned, or a tibble linking categories to stockGroups}

\item{stratAssign_fallback}{tibble with sWeek, stockGroup, stratum showing what stratum each sWeek
corresponds to for each stockGroup}

\item{stratAssign_comp}{tibble with sWeek, stratum showing what stratum each sWeek corresponds to
for composition estimation using the trap data}

\item{alpha_ci}{the alpha value for CI calculation}

\item{output_type}{one of "summary", "W_boot", or "full" dependign on the output type you want}

\item{split_HNC_fallbackone}{of "var1", "var2", or "both" indicating which variable(s)
should be used to link estimates of composition with fallback rates for "HNC" fish}

\item{split_W_fallbackone}{of "var1", "var2", or "both" indicating which variable(s)
should be used to link estimates of composition with fallback rates for "W" fish}
}
\description{
applies the estimated fallback rates to the estiamted ascension composition
to yield estimated escapement
}
