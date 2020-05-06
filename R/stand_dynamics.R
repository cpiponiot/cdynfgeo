#' Calculation of stand-level dynamics
#'
#' @param id A vector of individual trees (or stems) IDs
#' @param var A numerical vector of the same size as id, containing the
#'   measurements of the target variable
#' @param dvar Optional, a vector with variable change between 2 censuses.
#' @param year A numerical vector (same size as id), containing census years
#' @param group Optional; a vector (same size as id), containing a grouping factor.
#' @param site Optional; a vector (same size as id), containing the measurements' site
#'
#' @return A data.table (data.frame) with stocks, annual gain and loss of var
#'   for every census year and group level.
#'
#' @export
#'
#'
stand_dynamics = function(id,
                          var,
                          dvar = NULL,
                          year,
                          group = NULL,
                          plot = NULL,
                          site = NULL) {
  library(data.table)
  df = data.table(id, var, dvar, year, group, plot, site)
  if (is.null(group))
    df[, group := 1]
  if (is.null(plot))
    df[, plot := 1]
  if (is.null(site))
    df[, site := 1]

  df = subset(df,!is.na(var))

  all_years = unique(df[, c("year", "site")])
  dyears = all_years[, .(year = c(NA, sort(year)),
                         nextyear = c(sort(year), NA)), .(site)]
  df = merge(df, dyears, by = c("site", "year"))

  setorder(df, id, year)
  if (is.null(dvar))
    df[, dvar := c(diff(var) / diff(year), NA)]

  # last stem measurement: dvar = mortality
  df[, lastMeas := c(id[-1] != id[-nrow(df)], TRUE)]
  df[(lastMeas), dvar := var / (nextyear - year)]

  # gain and stocks
  dyn = df[, .(stock = sum(var), gain = sum(dvar[!lastMeas])), .(site, plot, group, year)]

  # mortality
  mort = df[, .(
    var = last(var),
    year = last(year),
    group = last(group)
  ), .(id, site, plot)]
  mort = merge(mort, dyears, by = c("site", "year"))
  mort[, mort := var / (nextyear - year)]
  df_loss = mort[, .(loss = sum(mort)), .(site, plot, group, year, nextyear)]
  dyn = merge(dyn, df_loss[,-"nextyear"],
              by = c("site", "plot", "group", "year"),
              all.x = TRUE)
  dyn = merge(dyn, dyears, by = c("site", "year"), all.x = TRUE)
  dyn[is.na(loss) & !is.na(nextyear), loss := 0]

  ## get all combinations of year * group per site
  combi = rbindlist(lapply(unique(df$site), function(s){
    expand.grid(site = s,
                plot = unique(df[site==s]$plot),
                group = unique(df[site==s]$group),
                year = unique(df[site==s]$year))
  }))
  dyn = merge(dyn, combi, by = c("site", "plot", "year", "group"), all = TRUE)
  # replace missing values by 0, except on last year (fluxes = NA)
  dyn[is.na(stock), stock := 0]
  dyn[is.na(gain), gain := 0]
  dyn[is.na(loss), loss := 0]
  dyn = merge(dyn[,-"nextyear"], dyears, by = c("site", "year"))
  dyn[is.na(nextyear), `:=`(gain=NA, loss=NA)]

  # remove unecessary columns
  dyn[, nextyear := NULL]
  if (is.null(group))
    dyn[, group := NULL]
  if (is.null(plot))
    dyn[, plot := NULL]
  if (is.null(site))
    dyn[, site := NULL]

  return(dyn)
}

# test

# load("temp_file.rda")
# data = fgeo_data[site=="scbi"]
#
# id = data$stemid
# var = rep(1, nrow(data))
# year = data$year
# group = as.numeric(cut(data$dbhc, 10))
# site = data$quadrat

# x=subset(fgeo_data, site=="bci"&quadrat=="4514" & dbhc < 100)
# id=x$stemid
# var=x$agb
# year=x$year
# group=x$size
# site=NULL
# dvar = NULL
# inst_flux = FALSE
