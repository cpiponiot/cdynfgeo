#' Calculation of stand-level dynamics (AGB, AWP, AWM, abundance, diameter
#' growth, stem mortality rate) used in the tree size analysis
#'
#' @param dbh A vector of numerical values, containing all dbh measurements.
#' @param agb A vector of numerical values, containing all agb measurements
#'   (same length as `dbh`)
#' @param year A vector of numerical values, containing census years
#'   corresponding to each measurement (same length as `dbh`).
#' @param stemid A vector of character (or numeric) values, containing unique
#'   stem identifiers corresponding to each measurement (same length as `dbh`).
#' @param hom A vector of numeric values, containing height of measurement, used
#'   to detect change in hom. Use `NA` when no HOM has been recorded.
#' @param wd A vector numeric values, containing wood density, used to
#'   (re-)calculate AGB
#' @param size A vector of character values corresponding to size bins
#' @param group Additional grouping, optional (default in `NULL`).
#' @param plot A vector of character (or numeric) values, containing subplot
#'   variable (quadrat or bigger) corresponding to each measurement (same length
#'   as `dbh`).
#' @param E A numerical value: environmental variable in Chave equation (2014)
#' @param a A numerical vector (same size as `dbh`): first parameter in
#'   Chojnacky equations (intercept). If tropical site, use `NA`.
#' @param b A numerical vector (same size as `dbh`): second parameter in
#'   Chojnacky equations (slope). If tropical site, use `NA`.
#' @param plot_area Area of subplots (in ha). Default is `1`.
#' @param Ddbh_range Vector of length 2: ange of acceptable DBH change (in cm). Default is c(-Inf, Inf) (ie no correction).
#'
#' @return A data.table (data.frame) with the following columns: `plot` and
#'   `year` and `size` as provided in the function inputs; `variable`: `N`
#'   (abundance, N/ha), `mrate` (mortality rate, %/yr), `Dgrowth` (dbh growth,
#'   cm/yr), `AGB` (aboveground biomass, Mg/ha), `AWP` (aboveground wood
#'   productivity, Mg/ha/yr), , `AWM` (aboveground wood mortality, Mg/ha/yr);
#'   `value` is the correponding value; `dT` is the time interval between the
#'   census and the following one (in years); `weight` is the weigth that will
#'   be used in the bootstrap function.
#'
#' @export
#'
Cdynamics = function(dbh,
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
                     Ddbh_range = c(-Inf, Inf)) {
  library(data.table)

  df = data.table(dbh, size, agb, year, stemid, hom, wd, plot, group, E, a, b)
  if (is.null(group))
    df$group <- size

  setorder(df, stemid, year)
  df[, dT := c(diff(year), NA)]
  # last stem measurement: dT = NA
  df[, lastMeas := c(stemid[-1] != stemid[-nrow(df)], TRUE)]
  df[(lastMeas), dT := NA]
  ## dbh change
  df[, Ddbh := c(diff(dbh), NA) / dT]
  # HOM change
  df[, dHOM := c(diff(hom), NA)]
  df[(lastMeas), dHOM := NA]
  df[, dHOM := (!is.na(dHOM) & dHOM != 0)]
  # when change in HOM: substitute with NA (substituted in the function below)
  # mean dbh change by site and size class
  ## DEBUG (if needed)
  # browser()
  # df[, .(
  #   pars = substitute_change(varD = Ddbh, hom_change = dHOM, value = "debug")
  # ), size]
  df[, `:=`(
    Ddbh = substitute_change(varD = Ddbh, hom_change = dHOM, cut = Ddbh_range),
    Dagb = substitute_change(
      D = dbh,
      varD = Ddbh,
      WD = wd,
      hom_change = dHOM,
      E = unique(E),
      a = a,
      b = b,
      value = "AGB",
      cut = Ddbh_range
    )
  ),
  size]

  ## abundance
  N_dyn = stand_dynamics(
    id = df$stemid,
    var = rep(1, nrow(df)),
    year = df$year,
    group = df$group,
    plot = df$plot
  )
  N_dyn[, N := stock]
  N_dyn[, mrate := loss]
  N_dyn = N_dyn[,-c("stock", "gain", "loss")]

  ## diameter growth
  D_dyn = stand_dynamics(
    id = df$stemid,
    var = df$dbh,
    dvar = df$Ddbh,
    year = df$year,
    group = df$group,
    plot = df$plot
  )
  D_dyn[, Dgrowth := gain]
  D_dyn = D_dyn[,-c("stock", "gain", "loss")]

  dyn = merge(N_dyn, D_dyn, by = c("year", "group", "plot"))

  ## biomass dynamics
  agb_dyn = stand_dynamics(
    id = df$stemid,
    var = df$agb,
    dvar = df$Dagb,
    year = df$year,
    group = df$group,
    plot = df$plot
  )
  colnames(agb_dyn)[4:6] = c("AGB", "AWP", "AWM")

  dyn = merge(dyn, agb_dyn, by = c("year", "group", "plot"))

  meas_vars = c("N", "mrate", "Dgrowth", "AGB", "AWP", "AWM")

  # dbh growth and stem mortality: standardized by number of stems
  n_vars = c("Dgrowth", "mrate")
  dyn[, wn := N]
  dyn[, n_vars] = dyn[, n_vars, with = FALSE] / dyn$wn
  # area_vars: weighted by area
  area_vars =  c("N", "AGB", "AWP", "AWM")
  dyn[, warea := plot_area]
  dyn[, area_vars] = dyn[, area_vars, with = FALSE] / dyn$warea
  dyn = melt(dyn, measure.vars = meas_vars, variable.factor = FALSE)

  ## add census interval dT to dyn
  dyear = data.table(year = sort(unique(df$year)))
  dyear = dyear[, dT := c(diff(year), NA)]
  # last census: dT = NA
  dyn = merge(dyn, dyear, by = "year")

  dyn = subset(dyn,!is.na(value))

  ## complete with zero value when there is no tree
  DF = expand.grid(
    group = unique(dyn$group),
    year = sort(unique(dyn$year), decreasing = TRUE)[-1],
    plot = unique(dyn$plot),
    variable = unique(dyn$variable)
  )
  dfdyn = merge(dyn,
                DF,
                by = c("group", "plot", "year", "variable"),
                all.y = TRUE)
  dfdyn[is.na(value), `:=`(value = 0, wn = 0, warea = plot_area)]

  # plot weight in bootstrap = standardizing variable (area or number of stems)
  dfdyn[, weight := warea]
  dfdyn[variable %in% n_vars, weight := wn]
  dfdyn = dfdyn[, -c("wn", "warea")]

  if (is.null(group)) {
    dfdyn[, `:=`(size = group, group = NULL)]
  }

  return(dfdyn)
}


## get symetrical distribution of dbh change by transforming with modulus function
## then calculate mean dbh change and backtransform


## modulus(Ddbh) not exactly normal: bimodal distribution (looks more like N(0,sd) + logN(mu, sd2))

# test
# x = fgeo_data[site=="bci"&dbhc>50]
# x[, dvar := c(diff(dbhc) / diff(year), NA)]
# # last stem measurement: dvar = NA
# x[, lastMeas := c(stemid[-1] != stemid[-nrow(x)], TRUE)]
# x[(lastMeas), dvar := NA]
# varD = x$dvar
# D = x$dbhc
# WD = x$wsg
# E = 1.5
# hom_change = x$dHOM
#

