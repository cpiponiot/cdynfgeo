#' ForestGEO-like data consolidation
#'
#' @param df A data.table obtained with the prepare_data() function.
#' @param taper_correction Logical value: should the taper correction from Cushman et al (2014) be applied? Default is `TRUE`.
#' @param stem_matching Logical value: should the stem matching script be applied? Default is `TRUE`.
#' @param add_missing Logical value: should missing stem measurements be interpolated? Default is `TRUE`.
#' @param correct_diam Logical value: should dbh errors (detected when dbh changes are
#'   outside the acceptable range) be corrected? Default is `TRUE`.
#' @param step_corr Logical value: should step dbh changes be corrected? Default is `TRUE`.
#' @param acc_decr numeric value: acceptable decrease between two censuses, as a proportion
#'   if relat_change is TRUE or in mm if relat_change is false. Default is -5 mm/yr
#' @param acc_incr numeric: acceptable annual increase, as a proportion (per year)
#'   if relat_change is TRUEor in mm/year if relat_change is false. Default is 35 mm/yr
#' @param relat_change logical value: should the decrease and increse values be relative
#'  (proportion of total dbh) or absolute values (in mm)? Default is `FALSE`.
#' @param species_path character string specifying the the directory in which the species tables are. Default is the current directory.
#'
#' @return A data.table (data.frame) with all relevant variables.
#'
#' @export
#'

consolidate_data <- function(df,
                             taper_correction = TRUE,
                             stem_matching = TRUE,
                             add_missing = TRUE,
                             correct_diam = TRUE,
                             step_corr = TRUE,
                             dcor_min = 100,
                             acc_decr = -5,
                             acc_incr = 50,
                             relat_change = FALSE,
                             species_path = getwd()) {

  library(data.table)
  #### remove unecessary rows ####
  t0 = Sys.time()
  df$status = c("D", "A")[1 + grepl("alive|broken", df$dfstatus)]
  print("Changed status to A (alive) or D (dead).")
  print(round(Sys.time() - t0, 1))
  t0 = Sys.time()

  df = subset(df, !is.na(dbh) & dbh >= 10 & status == "A")
  print("Removed non-observations (no dbh, dead stems, stems < 10 mm).")
  print(round(Sys.time() - t0, 1))
  t0 = Sys.time()

  #### make unique tree info table ####
  ## keep the most recent info
  df$treeid = as.character(df$treeid)
  treeinfo = df[order(year),
                .(quadrat = last(quadrat[!is.na(quadrat)]),
                  gx = last(gx[!is.na(gx)]),
                  gy = last(gy[!is.na(gy)]),
                  sp = last(sp[!is.na(sp)])),
                by = treeid]
  df = merge(df[,-c("quadrat", "gx", "gy", "sp")], treeinfo, by = "treeid")

  #### add accepted species name ####
  if (!("species_clean.rda" %in% list.files(species_path))) {
    # get clean version of species table (no typo nor duplicate species code)
    species = clean_spptab(species_path)
    save(species, file = paste0(species_path, "species_clean.rda"))
  } else
    load(paste0(species_path, "species_clean.rda"))
  df = merge(df, species, by = c("sp", "site"), all.x = TRUE)
  # df[, genus := tstrsplit(name, " ")[[1]]]
  # df[, species := tstrsplit(name, " ")[[2]]]
  print("Added species name.")
  print(round(Sys.time() - t0, 1))
  t0 = Sys.time()

  #### taper correction ####
  df[, dbhc := dbh] # dbhc: corrected dbh
  if (taper_correction) {
    df[, dbhtc := taper(dbh, hom)]
    df[, dbhc := dbhtc]
    print("Applied taper correction.")
    print(round(Sys.time() - t0, 1))
    t0 = Sys.time()
    # add uncertainty on dbh measurement to propagate later
  }

  #### match stems - creates new stem ID ####
  if (stem_matching) {
    df[, stemid := paste(treeid, stemtag, sep = "_")]
    df[!is.na(stemid), stemYearID := paste(treeid, stemid, year, sep = "_")]
    # determine which stems need matching: no stem tag or duplicated stem tag
    need_stemmatch = unique(df[is.na(stemtag) |
                                 duplicated(stemYearID), treeid])
    # apply stem matching and merge with df based on rowid
    df[, rowid := (1:nrow(df))]
    newstem = df[treeid %in% need_stemmatch,
                 .(
                   stemid2 = match_stems(dbhc, year, acc_decr, acc_incr, relat_change),
                   rowid = rowid
                 ),
                 .(treeid)]
    df = merge(df,
               newstem,
               by = c("treeid", "rowid"),
               all = TRUE)
    df[treeid %in% need_stemmatch, stemid := paste(treeid, stemid2, sep =
                                                     "_")]
    print("Matched stems in multistem individuals (in column stemid).")
    print(round(Sys.time() - t0, 1))
    t0 = Sys.time()

    ## broken below: new stem id (add 'b' at the end for all subsequent censuses)
    setorder(df, stemid, year)
    # time interval between censuses (to detect overgrown resprouts)
    df[, dt := c(NA, diff(year))]
    df[c(NA, stemid[-length(stemid)])!= c(NA, stemid[-1]), dt := NA]
    # real resprouts:
    broken = df[codes == "R" & dbhc <= acc_incr*dt, c("stemid", "year")]
    several_breaks = TRUE
    if (several_breaks & nrow(broken) > 0) {
      broken = broken[, .(nbreak = paste0("break", 1:length(year)), year = sort(year)), .(stemid)]
      broken = dcast(broken, stemid ~ nbreak, value.var = "year")
      df = merge(df, broken, by = c("stemid"), all = TRUE)
      df[!is.na(break1) & year >= break1, stemid := paste0(stemid, "b")]
      if (!is.null(df$break2))
        df[!is.na(break2) & year >= break2, stemid := paste0(stemid, "b")]
      if (!is.null(df$break3))
        df[!is.na(break3) & year >= break3, stemid := paste0(stemid, "b")]
    } else if (nrow(broken) > 0) {
      broken = broken[, .(ybreak = min(year)), .(stemid)]
      df = merge(df, broken, by = c("stemid"), all = TRUE)
      df[!is.na(ybreak) & year >= ybreak, stemid := paste0(stemid, "b")]
    }
  }

  ### add growth statistics ####
  dfgrowth = growth_stats(
    dbh = df$dbhc,
    year = df$year,
    stemid = df$stemid,
    site = df$site,
    name = df$name,
    relat_change = relat_change
  )
  # expected growth rate (lower taxonimic level available)
  df = merge(df, dfgrowth, by = c("site", "name"), all.x = TRUE)

  #### add missing tree dbh ####
  if (add_missing) {
    t0 = Sys.time()
    # add missing measurements (all censuses between recruitment and death)
    df_census = df[, .(census = min(census):max(census)), .(stemid, site)]
    df_census[, status := "A"]
    # add dbh values (dbh = NA: stem missing)
    df_census = merge(df_census,
                      df[, c("stemid", "census", "dbhc")],
                      by = c("stemid", "census"),
                      all = TRUE)
    # complete with census year
    census_year = unique(df[, c("site", "census", "year")])
    df_census = merge(df_census, census_year , by = c("site", "census"))

    # interpolate missing dbhs
    setorder(df_census, year)
    x = df_census[, .(dbhc = replace_missing(dbhc, year, dfgrowth), year), .(stemid)]
    df = merge(df[,-"dbhc"], x, by = c("stemid", "year"), all = TRUE)
    print("Interpolated DBH for missing trees.")
    print(round(Sys.time() - t0, 1))
  }

  #### correct abnormal dbh changes ####
  if (correct_diam) {
    t0 = Sys.time()
    df2 = df[, .(dbhcc = correct_dbh(dbh = dbhc,
                                     year = year,
                                     Gexp = unique(Gexp[!is.na(Gexp)]),
                                     codes = codes,
                                     acc_decr = acc_decr,
                                     acc_incr = acc_incr,
                                     relat_change = relat_change,
                                     dcor_min = dcor_min,
                                     step_corr = step_corr),
                 year),
             .(stemid)]
    df = merge(df, df2, by = c("stemid", "year"))
    df[, dbhc := dbhcc]
    df[dbhc < 10, dbhc := NA]
    print("Corrected excessive dbh changes.")
    print(round(Sys.time() - t0, 1))
    t0 = Sys.time()
  }

  # #### merge tree and census info ####
  df_tree = df[, c("site", "quadrat", "treeid", "gx", "gy",
                   "sp", "name", "taxo_group")]
  df_tree = subset(df_tree, !is.na(sp))
  df_tree = unique(df_tree)

  cols = c("treeid", "stemid", "year", "hom", "dbh", "dbhc", "dbhtc", "dbhcc")
  df2 = df[, which(colnames(df) %in% cols), with = FALSE]
  df = merge(df_tree, df2, by = "treeid")

  return(df)
}


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



# growth_stats = function(dbh,
#                         year,
#                         stemid,
#                         site,
#                         name,
#                         relat_change = FALSE) {
#   data = data.table(site, stemid, name, dbh, year)
#   data[, genus := tstrsplit(name, " ")[[1]]]
#   data[, species := tstrsplit(name, " ")[[2]]]
#
#   # remove individuals with only one measurement (no diff(dbh))
#   multi = unique(stemid[duplicated(stemid)])
#   setorder(data, stemid, year)
#   dataG = data[stemid %in% multi]
#   if (relat_change) {
#     dataG = dataG[, .(dG = (dbh[-1] / dbh[-length(dbh)]) ^ (1 / diff(year)) -
#                         1,
#                       year = year[-1]),
#                   .(stemid, name, genus, species, site)]
#   } else {
#     dataG = dataG[, .(dG = diff(dbh) / diff(year), year = year[-1]),
#                   .(stemid, name, genus, species, site)]
#   }
#
#   dfnames = unique(data[,c("site", "name", "genus")])
#
#   # species growth rate (>20 individual measurements)
#   dfspecies = dataG[!is.na(species), .(
#     Gexp = median(dG), ## use median to minimise impact of measurement errors
#     sdGexp = sd(dG),
#     N = length(dG),
#     level = "species"
#   ),
#   .(site, genus, name)]
#   dfspecies = subset(dfspecies, N >= 20)[, -"N"]
#
#   # genus growth rate (>20 individual measurements)
#   dfgenus = dataG[!is.na(genus), .(
#     Gexp = median(dG),
#     sdGexp = sd(dG),
#     N = length(dG),
#     level = "genus"
#   ),
#   .(genus)]
#   dfgenus = subset(dfgenus, N >= 20)[, -"N"]
#   dfgenus = merge(dfgenus, dfnames, by = "genus")
#
#   dfsite = dataG[, .(Gexp = median(dG),
#                      sdGexp = sd(dG),
#                      level = "site"), .(site)]
#   dfsite = merge(dfsite, dfnames, by = "site")
#
#   dfgrowth = merge(dfnames, rbind(dfspecies, dfgenus, dfsite), all = TRUE)
#   dfgrowth[is.na(Gexp)]$Gexp = mean(dataG$dG)
#   dfgrowth[is.na(sdGexp)]$sdGexp = sd(dataG$dG)
#   dfgrowth = dfgrowth[!duplicated(dfgrowth[, c("site", "name")])]
#
#   return(dfgrowth)
# }
