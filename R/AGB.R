#' Calculate aboveground biomass for ForestGEo sites
#'
#' @param dbh A vector of numerical values, containing all dbh measurements.
#' @param wd A vector numeric values, containing wood density (same length as
#'   `dbh`).
#' @param E A numerical value: environmental variable in Chave equation (2014)
#' @param a A numerical vector (same size as `dbh`): first parameter in
#'   Chojnacky equations (intercept). If tropical site, use `NA`.
#' @param b A numerical vector (same size as `dbh`): second parameter in
#'   Chojnacky equations (slope). If tropical site, use `NA`.
#' @param carbon Should the function return values of aboveground carbon
#'   (instead of biomass)? Uses the conversion factor from Martin et al 2012.
#'   Default is `FALSE`.
#'
#' @return A vector of aboveground biomass (or carbon if `carbon = TRUE`).
#'
#' @export
#'
AGB = function(dbh, wd, E, a, b, carbon = FALSE) {
  ### tropical sites: Chave equation
  agb = exp(
    -2.023977 - 0.89563505 * E + 0.92023559 *
      log(wd) + 2.79495823 * log(dbh) - 0.04606298 * (log(dbh) ^
                                                        2)
  )

  ## temperate sites: Chojnacky equations
  agb[!is.na(a) &
        !is.na(b)] = exp(a + b * log(dbh))[!is.na(a) & !is.na(b)]

  if (carbon) {
    agb[!is.na(a) &
          !is.na(b)] = agb[!is.na(a) &
                             !is.na(b)] * 0.49 ## include if possible distrinction conifer/angiosperm
    agb[is.na(a) & is.na(b)] = agb[is.na(a) & is.na(b)] * 0.471
  }

  return(agb)

}
