#' Substitute DBH changes (cm/yr) outside of a pre-set acceptable range of DBH
#' changes by their expected values based on other observed values inside the
#' pre-set interval. Can also substitute AGB changes by intergrating the AGB
#' allometry over the distribution of expected DBH changes. This function allows
#' to have unbiased stand-level estimates of DBH and AGB growth while removing
#' the most obvious measurement errors.
#'
#' @param varD Vector of changes in DBH.
#' @param hom_change Vector of changes in height of measurement.
#' @param cut Acceptable range of DBH changes. Default is c(-0.5, 5)
#' @param lambda Parameter of the modulus function used to normalize DBH changes
#'   (see Condit et al., 2017). Default is 0.5
#' @param value What change should be computed, either "D" (DBH) or "AGB".
#'   Default is "D".
#' @param D Must provide a value when value = "AGB". Default is NULL.
#' @param WD Vector of wood densities. Must provide a value when value = "AGB"
#'   (can be "NA" for temperate sites). Default is NULL.
#' @param E Vector of environmental parameter used in Chave et al pantropical
#'   allometric equation. Must provide a value when value = "AGB" (NA for
#'   temperate sites). Default is NULL.
#' @param a Vector of first parameter used in Chojnacky allometric equations.
#'   Must provide a value when value = "AGB" (NA for tropical sites). Default is
#'   NULL.
#' @param b Vector of second parameter used in Chojnacky allometric equations.
#'   Must provide a value when value = "AGB" (NA for tropical sites). Default is
#'   NULL.
#'
#' @return A vector of corrected changes in DBH (if value = "D) or AGB (if value
#'   = "AGB).
#'
#' @export
#'
substitute_change = function(varD,
                             hom_change,
                             cut = c(-0.5, 5),
                             lambda = 0.5,
                             value = "D",
                             D = NULL,
                             WD = NULL,
                             E = NULL,
                             a = NULL,
                             b = NULL) {

  if (value == "AGB" & (is.null(D) | is.null(WD) | is.null(E) | is.null(a) | is.null(b)))
    stop("You have to provide D, WD, E, a and b in order to correct AGB changes.")

  keep_values = which(!hom_change &
                        varD > cut[1] & varD < cut[2] & !is.na(varD))
  change_values = which((hom_change |
                           varD <= cut[1] | varD >= cut[2]) & !is.na(varD))
  transf_values = modulus(varD[keep_values], lambda)
  mu = mean(transf_values)
  sigma = sd(transf_values)

  if (value == "debug")
    return(c(mu, sigma))

  if (value == "D") {
    diffModD = function(x)
      modulus(x, 1 / lambda) * dnorm(x, mu, sigma)
    varD[change_values] = integrate(diffModD, -Inf, Inf)$value
    return(varD)
  }

  if (value == "AGB") {
    dAGB = AGB(D + varD, WD, E, a , b) - AGB(D, WD, E, a, b)
    if (length(change_values) == 1) {
      dAGB[change_values] = ExpDiffAGB(
        d = D[change_values],
        wd = WD[change_values],
        E = E,
        a = a[change_values],
        b = b[change_values],
        lambda = lambda,
        mu = mu,
        sigma = sigma
      )
    } else if (length(change_values) > 1) {
      dAGB[change_values] = apply(cbind(D, WD, a, b)[change_values, ], 1, function(x) {
        ExpDiffAGB(
          d = x[1],
          wd = x[2],
          E = E,
          a = x[3],
          b = x[4],
          lambda = lambda,
          mu = mu,
          sigma = sigma
        )
      })
    }
    return(dAGB/1000)
  }
}

modulus = function(d, lambda = 0.4) {
  return(sign(d) * abs(d) ^ lambda)
}

ExpDiffAGB = function(d, wd, E, a, b, lambda, mu, sigma) {
  # AGB > 0 => d + modulus(x, 1/lambda) > 0 => x > -(d^lambda)
  minVar = -d ^ lambda
  pdfdAGB = function(x) {
    dAGB = AGB(d + modulus(x, 1 / lambda), wd, E, a, b) - AGB(d, wd, E, a, b)
    dens = truncnorm::dtruncnorm(x, mean = mu, sd = sigma, a = minVar)
    return(dAGB * dens)
  }
  return(integrate(pdfdAGB, lower = minVar, upper = Inf)$value)
}
