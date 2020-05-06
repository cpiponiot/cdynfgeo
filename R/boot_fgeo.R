#' Bootstrap values to get site-level estimates and 95% confidence intervals
#'
#' @param value Data to be bootstraped, as a numerical vector.
#' @param weight Weight given to observations, as a numerical vector of same
#'   length as `value`.
#' @param subplot Subplots observations belong to.
#' @param variable Variables observations belong to.
#' @param group Grouping factor (for example size).
#' @param group2 Additional grouping factor (for example site).
#' @param year Year observations belong to. Default is `NULL`.
#' @param dT Time interval (years) between the current census (as in `year`) and
#'   the next one.
#' @param kohyama Name of variables to apply Kohyama correction (2019), should
#'   be a vector of length 3 with name of stock variable, influx variable, and
#'   outflux variable (in that order). Default is `NULL` (i.e. correction not
#'   applied).
#' @param nrep Default is 1000.
#' @param para Default is `FALSE`.
#'
#' @return A data frame with relevant variables.
#'
#' @export
#'
boot_fgeo = function(value,
                     weight,
                     subplot,
                     variable,
                     group,
                     group2,
                     year = NULL,
                     dT = NULL,
                     kohyama = NULL,
                     nrep = 1000,
                     para = FALSE) {
  data = data.table(value, weight, group, subplot, year, variable, group2, dT)
  # parallelized (or not) bootstrapping
  if (para) {
    library(parallel)
    cl <- makeCluster(detectCores() - 1)
    Ldf = data.table::split(data, by = "group2", drop = TRUE)
    clusterEvalQ(cl, library(data.table)) # packages
    clusterExport(cl, varlist = c("Ldf", "nrep", "bootstrap", "kohyama_correction", "kohyama")) # data
    listboot <- parLapply(cl, Ldf, function(x) bootstrap(x, nrep, kohyama))
    stopCluster(cl)
  } else {
    Ldf = split(data, by = "group2", drop = TRUE)
    listboot <- lapply(Ldf, function(x)
      bootstrap(x, nrep, kohyama))
  }
  # get group2 and variable from list names
  for (i in 1:length(listboot))
    listboot[[i]]$group2 = names(listboot)[i]
  # results
  DF = rbindlist(listboot)
  return(DF)
}

bootstrap = function(df, n, var_corr) {
  # split into subplots
  Lobs = split(df, by = "subplot")
  # sample observations (with replacement)
  samples = sample(length(Lobs), n * length(Lobs), replace = TRUE)
  sample_reps = split(samples, rep(1:n, each = length(Lobs)))
  # create new resampled data frames and calculate weigthed mean
  samplist = lapply(sample_reps, function(i) {
    dt = rbindlist(Lobs[i])
    if (length(var_corr) > 0) {
      if (sum(var_corr %in% dt$variable) != 3 | length(df$dT) == 0)
        stop("You should provide 3 existing variables to apply the Kohyama correction, and census interval dT.")
      dt = kohyama_correction(dt, var_corr)
    }
    # get weighted mean
    xw = dt$value * dt$weight
    group_var = paste(dt$group, dt$variable)
    resu = tapply(xw, group_var, sum) / tapply(dt$weight, group_var, sum)
    new_df = data.frame(
      x = resu,
      group = tstrsplit(names(resu), " ")[[1]],
      variable = tstrsplit(names(resu), " ")[[2]]
    )
    return(new_df)
  })
  sampled = rbindlist(samplist)
  ci = sampled[, .(lwr = quantile(x, 0.025, na.rm = TRUE),
                   upr = quantile(x, 0.975, na.rm = TRUE)),
               .(group, variable)]
  xwtot = tapply(df$value * df$weight, paste(df$group, df$variable), sum)
  wtot = tapply(df$weight, paste(df$group, df$variable), sum)
  ci$all = xwtot / wtot
  return(ci)
}

# instanteneous biomass (or other) fluxes as recommended by Kohyama 2019 (eq 1-2 in Table 1)
kohyama_correction = function(dt, vars){
  dt_sub = dt[variable %in% vars, .(value = sum(value*weight)/sum(weight), weight = sum(weight)),
              .(dT, year, variable, group, group2)]
  dt_new = dcast(dt_sub, year + dT + group2 + weight + group ~ variable)
  dt_new$B0 = dt_new[, vars[1], with = FALSE]
  dt_new$gain = dt_new[, vars[2], with = FALSE]
  dt_new$loss = dt_new[, vars[3], with = FALSE]
  dt_new[, BS0 := B0 - loss*dT]
  dt_new[, BT := B0 + (gain-loss)*dT]
  dt_new[, vars[2] := (log(BT/BS0)*(BT-B0))/(dT*log(BT/B0))]
  dt_new[, vars[3] := (log(B0/BS0)*(BT-B0))/(dT*log(BT/B0))]

  dt_new = melt(dt_new, id.vars = c("group2", "weight", "group"), measure.vars = vars)
  dt_corr = rbind(dt[!variable %in% vars, colnames(dt_new), with = FALSE], dt_new)
  dt_corr = subset(dt_corr, !is.na(value))
  return(dt_corr)
}

