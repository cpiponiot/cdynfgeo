boot_fgeo = function(value,
                     weight,
                     group,
                     subplot,
                     variable,
                     site,
                     year = NULL, 
                     dT = NULL, 
                     kohyama = NULL,
                     nrep = 1000,
                     para = FALSE) {
  data = data.table(value, weight, group, subplot, year, variable, site, dT)
  # parallelized (or not) bootstrapping
  if (para) {
    library(parallel)
    cl <- makeCluster(detectCores() - 1)
    Ldf = split(data, by = "site", drop = TRUE)
    clusterEvalQ(cl, library(data.table)) # packages
    clusterExport(cl, varlist = c("Ldf", "nrep", "bootstrap", "kohyama_correction", "kohyama")) # data
    listboot <- parLapply(cl, Ldf, function(x) bootstrap(x, nrep, kohyama))
    stopCluster(cl)
  } else {
    Ldf = split(data, by = "site", drop = TRUE)
    listboot <- lapply(Ldf, function(x)
      bootstrap(x, nrep, kohyama))
  }
  # get site and variable from list names
  for (i in 1:length(listboot))
    listboot[[i]]$site = names(listboot)[i]
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
