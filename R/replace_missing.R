replace_missing <- function(dbh, year, Gexp) {
  
  # TODO check that there are at least 2 measurements
  miss <- which(is.na(dbh))
  nomiss <- which(!is.na(dbh))
  
  if (length(miss) == 0){
    return(dbh)
  } else if (length(nomiss)==1) {
    # Gexp: expected growth rate (based on species- 
    # or genus- or site- specific growth rates)
    dbh[miss] = Gexp*(year[miss]-year[nomiss]) + dbh[nomiss]
    return(dbh)
    
  } else {
    # for each missing value, select 2 closest measurements 
    # from which missing dbh value will be inferred
    inf_years = sapply(miss, function (i) select_2years(i, nomiss))
    
    xs = matrix(year[inf_years], nrow=2)
    ys = matrix(dbh[inf_years], nrow=2)
    
    dbh[miss] = interpolate(year[miss], xs, ys) 
  }
  return(dbh)
}

select_2years = function(i, nomiss) {
  if (i < min(nomiss)){ 
    return(nomiss[1:2])
  } else if (i > max(nomiss)) {
    last = length(nomiss)
    return(nomiss[(last-1):last])
  } else {
    return(c(max(nomiss[nomiss<i]),  min(nomiss[nomiss>i])))
  }
}

interpolate = function(x, xs, ys) {
  slope = diff(ys)/diff(xs)
  intercept = ys[1,] - slope*xs[1,]
  return(slope * x + intercept)
}

