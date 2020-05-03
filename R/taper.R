taper <- function(dbh, hom, wsg = NULL) {
  hom[is.na(hom)] = 1.3
  # dbh <- round(dbh * exp(0.0247 * (hom - 1.3)), 1)
  # new Cushman multisite equation
  # DAB: diameter above buttress / same as dbh?
  if (is.null(wsg)) {
    b1 <- exp(-0.825 - 0.458 * log(dbh / 10) - 0.663 * log(hom))
  } else {
    b1 <- exp(-0.769 - 0.546 * log(dbh / 10) - 0.630 * log(hom) - 0.436 * log(wsg))
  }
  # b1 = -0.825 + -0.458*log(DAB) -0.663*log(HOM)
  dbhc <- round(dbh * exp(b1 * (hom - 1.3)), 1)
  # add uncertainty -> uncertainty on dbh -> error propagation in AGBmontecarlo
  # change hom to 1.3?
  return(dbhc)
}
