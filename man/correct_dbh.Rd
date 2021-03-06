% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct_dbh.R
\name{correct_dbh}
\alias{correct_dbh}
\title{DBH correction in ForestGEO data}
\usage{
correct_dbh(
  dbh,
  year,
  Gexp = 5,
  codes = NA,
  acc_decr = -5,
  acc_incr = 35,
  relat_change = FALSE,
  step_corr = TRUE,
  dcor_min = 10
)
}
\arguments{
\item{dbh}{a numeric vector containing the dbh measurements of a tree or stem}

\item{year}{a numeric vector containing the census years of a tree or stem}

\item{Gexp}{expected growth rate}

\item{codes}{optional, a character vector containing the field codes of a tree or stem}

\item{acc_decr}{numeric: acceptable decrease between two censuses, as a proportion
if relat_change is TRUE or in mm if relat_change is false}

\item{acc_incr}{numeric: acceptable annual increase, as a proportion (per year)
if relat_change is TRUEor in mm/year if relat_change is false}

\item{relat_change}{logical: should the decrease and increse values be relative
(proportion of total dbh) or absolute values (in mm)?}

\item{step_corr}{logical: should step dbh changes be corrected?}

\item{dcor_min}{numeric: minimum dbh (in mm) of corrected stems}
}
\value{
A numeric vector with corrected dbhs.
}
\description{
DBH correction in ForestGEO data
}
