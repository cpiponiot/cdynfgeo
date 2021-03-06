% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cdynamics.R
\name{Cdynamics}
\alias{Cdynamics}
\title{Calculation of stand-level dynamics (AGB, AWP, AWM, abundance, diameter
growth, stem mortality rate) used in the tree size analysis}
\usage{
Cdynamics(
  dbh,
  agb,
  year,
  stemid,
  hom,
  wd,
  size,
  group = NULL,
  plot,
  E,
  a,
  b,
  plot_area = 1,
  Ddbh_range = c(-Inf, Inf)
)
}
\arguments{
\item{dbh}{A vector of numerical values, containing all dbh measurements.}

\item{agb}{A vector of numerical values, containing all agb measurements
(same length as \code{dbh})}

\item{year}{A vector of numerical values, containing census years
corresponding to each measurement (same length as \code{dbh}).}

\item{stemid}{A vector of character (or numeric) values, containing unique
stem identifiers corresponding to each measurement (same length as \code{dbh}).}

\item{hom}{A vector of numeric values, containing height of measurement, used
to detect change in hom. Use \code{NA} when no HOM has been recorded.}

\item{wd}{A vector numeric values, containing wood density, used to
(re-)calculate AGB}

\item{size}{A vector of character values corresponding to size bins}

\item{group}{Additional grouping, optional (default in \code{NULL}).}

\item{plot}{A vector of character (or numeric) values, containing subplot
variable (quadrat or bigger) corresponding to each measurement (same length
as \code{dbh}).}

\item{E}{A numerical value: environmental variable in Chave equation (2014)}

\item{a}{A numerical vector (same size as \code{dbh}): first parameter in
Chojnacky equations (intercept). If tropical site, use \code{NA}.}

\item{b}{A numerical vector (same size as \code{dbh}): second parameter in
Chojnacky equations (slope). If tropical site, use \code{NA}.}

\item{plot_area}{Area of subplots (in ha). Default is \code{1}.}

\item{Ddbh_range}{Vector of length 2: ange of acceptable DBH change (in cm). Default is c(-Inf, Inf) (ie no correction).}
}
\value{
A data.table (data.frame) with the following columns: \code{plot} and
\code{year} and \code{size} as provided in the function inputs; \code{variable}: \code{N}
(abundance, N/ha), \code{mrate} (mortality rate, \%/yr), \code{Dgrowth} (dbh growth,
cm/yr), \code{AGB} (aboveground biomass, Mg/ha), \code{AWP} (aboveground wood
productivity, Mg/ha/yr), , \code{AWM} (aboveground wood mortality, Mg/ha/yr);
\code{value} is the correponding value; \code{dT} is the time interval between the
census and the following one (in years); \code{weight} is the weigth that will
be used in the bootstrap function.
}
\description{
Calculation of stand-level dynamics (AGB, AWP, AWM, abundance, diameter
growth, stem mortality rate) used in the tree size analysis
}
