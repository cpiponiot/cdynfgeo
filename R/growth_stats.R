growth_stats = function(dbh,
                        year,
                        stemid,
                        site,
                        name,
                        relat_change = FALSE) {
  data = data.table(site, stemid, name, dbh, year)
  data[, genus := tstrsplit(name, " ")[[1]]]
  data[, species := tstrsplit(name, " ")[[2]]]
  
  # remove individuals with only one measurement (no diff(dbh))
  multi = unique(stemid[duplicated(stemid)])
  setorder(data, stemid, year)
  dataG = data[stemid %in% multi]
  if (relat_change) {
    dataG = dataG[, .(dG = (dbh[-1] / dbh[-length(dbh)]) ^ (1 / diff(year)) -
                        1,
                      year = year[-1]),
                  .(stemid, name, genus, species, site)]
  } else {
    dataG = dataG[, .(dG = diff(dbh) / diff(year), year = year[-1]),
                  .(stemid, name, genus, species, site)]
  }
  
  dfnames = unique(data[,c("site", "name", "genus")])
  
  # species growth rate (>20 individual measurements)
  dfspecies = dataG[!is.na(species), .(
    Gexp = median(dG), ## use median to minimise impact of measurement errors
    sdGexp = sd(dG),
    N = length(dG),
    level = "species"
  ),
  .(site, genus, name)]
  dfspecies = subset(dfspecies, N >= 20)[, -"N"]
  
  # genus growth rate (>20 individual measurements)
  dfgenus = dataG[!is.na(genus), .(
    Gexp = median(dG),
    sdGexp = sd(dG),
    N = length(dG),
    level = "genus"
  ),
  .(genus)]
  dfgenus = subset(dfgenus, N >= 20)[, -"N"]
  dfgenus = merge(dfgenus, dfnames, by = "genus")
  
  dfsite = dataG[, .(Gexp = median(dG),
                     sdGexp = sd(dG),
                     level = "site"), .(site)]
  dfsite = merge(dfsite, dfnames, by = "site")
  
  dfgrowth = merge(dfnames, rbind(dfspecies, dfgenus, dfsite), all = TRUE)
  dfgrowth[is.na(Gexp)]$Gexp = mean(dataG$dG)
  dfgrowth[is.na(sdGexp)]$sdGexp = sd(dataG$dG)
  dfgrowth = dfgrowth[!duplicated(dfgrowth[, c("site", "name")])]
  
  return(dfgrowth)
}