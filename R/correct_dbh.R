#' DBH correction in ForestGEO data
#'
#' @param dbh a numeric vector containing the dbh measurements of a tree or stem
#' @param year a numeric vector containing the census years of a tree or stem
#' @param Gexp expected growth rate
#' @param codes optional, a character vector containing the field codes of a tree or stem
#' @param acc_decr numeric: acceptable decrease between two censuses, as a proportion
#'   if relat_change is TRUE or in mm if relat_change is false
#' @param acc_incr numeric: acceptable annual increase, as a proportion (per year)
#'   if relat_change is TRUEor in mm/year if relat_change is false
#' @param relat_change logical: should the decrease and increse values be relative
#'  (proportion of total dbh) or absolute values (in mm)?
#' @param step_corr logical: should step dbh changes be corrected?
#' @param dcor_min numeric: minimum dbh (in mm) of corrected stems
#'
#' @return A numeric vector with corrected dbhs.
#'
#' @export

# test
# id = "bci_741_1"
# dbh = fgeo_data[stemid == id]$dbhtc
# year = fgeo_data[stemid == id]$year
# codes = fgeo_data[stemid == id]$codes
# # codes = NA
# Gexp = 2

correct_dbh <- function(dbh,
                        year,
                        Gexp = 5,
                        codes = NA,
                        acc_decr = -5,
                        acc_incr = 35,
                        relat_change = FALSE,
                        step_corr = TRUE,
                        dcor_min = 10) {
  # keep order of censuses (to return ordered dbhs)
  yorder = order(year)
  # order censuses
  dbhc = dbh[yorder]
  year = sort(year)
  # annual growth
  noNA = which(!is.na(dbhc))
  dG = anGrowth(dbhc[noNA], year[noNA], relat_change)
  # problematic growth rates
  prob = which(dG > acc_incr | dG < acc_decr)
  
  #### (1) no correction needed ####
  if (length(dG) == 0 |
      all(dG < acc_incr & dG > acc_decr) |
      all(dbhc[noNA[c(prob, prob + 1)]] < dcor_min)) {
    return(dbhc[order(yorder)])
    
  } else {
    
    #### (2) there are only two non missing measurements ####
    if (length(noNA) == 2) {
      if (! step_corr | 
          (dG < 0 &
           length(codes) > 1 &
           !is.na(codes[noNA[2]]) &
           codes[noNA[2]] == "R" &
           dbhc[noNA[2]] < acc_incr * diff(year))){
        ## if step_corr = FALSE: do nothing
        ## if the field code is "R" and dbh not too large: 
        ## it is a resprout and should not be corrected
        return(dbhc[order(yorder)])
      } else {
        if (dbhc[noNA[1]] > 100 & dG < 0) {
          # if the problem is the dbh decrease of a large tree
          # it is likely caused by a change in hom
          # then correct the more recent measurements
          dbhc[noNA[2]] = dbhc[noNA[1]] + Gexp * diff(year[noNA])
        } else {
          # else, trust the more recent measurement
          dbhc[noNA[1]] = dbhc[noNA[2]] - Gexp * diff(year[noNA])
        }
        return(dbhc[order(yorder)])
      }
      
    } else {
      #### (3) at least 3 dbhs and 1 error ####
      ## sequentially correct dbh errors, starting with the highest difference in dbh
      toCorr = which((dG > acc_incr | dG < acc_decr) & !is.na(dG))
      toCorr = toCorr[order(abs(dG[toCorr]), decreasing = TRUE)]
      
      for (corr in toCorr) {
        #### (3.a) a decrease caused by a broken and resprouting stem ####
        ## do not correct ##
        ## else: ##
        if (!(dG[corr] < 0 &
              !all(is.na(codes)) &
              !is.na(codes[noNA[corr + 1]]) &
              codes[noNA[corr + 1]] == "R" &
              dbhc[noNA[corr + 1]] < acc_incr * diff(year)[noNA[corr]]) &
            # check that this error has not been corrected yet
            (dG[corr] > acc_incr | dG[corr] < acc_decr)) {
          #### (3.b) one outlier ####
          # look for outliers
          # select surrounding values, correcting the dbh before or after
          # the abnormal dbh change
          diffBefore = noNA[c(corr - 1, corr + 1)]
          diffAfter = noNA[c(corr, corr + 2)]
          dGB = anGrowth(dbhc[diffBefore], year[diffBefore], relat_change)
          if (length(dGB) == 0)
            dGB = NA
          dGA = anGrowth(dbhc[diffAfter], year[diffAfter], relat_change)
          # choose the dbh change that is closer to expected change Gexp
          choice = which.min(c(abs(dGB - Gexp), abs(dGA - Gexp)))
          chdG = c(dGB, dGA)[choice]
          
          if (length(chdG) == 1 &
              ((chdG < acc_incr &
                chdG > acc_decr) | chdG / dG[corr] < 0.2)) {
            # remove erroneous dbh, and interpolate new value
            dbhc[noNA[corr + choice - 1]] = NA
            dbhc[noNA] = replace_missing(dbhc[noNA], year[noNA], Gexp)
            
          } else if (step_corr) {
            #### (3.c) step change in dbh ####
            # if the problem is a setp dbh change of a large tree
            # it is likely caused by a change in hom
            # then we correct the more recent measurements;
            # else we trust the set of measurements with more values, or if
            # if they are the same size, then we trust the more recent one
            if (corr > length(noNA) / 2 |
                dbhc[noNA[corr]] > 100) {
              dbhc[noNA[(corr + 1):length(noNA)]] <-
                dbhc[noNA[(corr + 1):length(noNA)]] -
                diff(dbhc[noNA[corr + 0:1]]) + Gexp * diff(year[noNA[corr +
                                                                       0:1]])
            } else {
              dbhc[noNA[1:corr]] <- dbhc[noNA[1:corr]] +
                diff(dbhc[noNA[corr + 0:1]]) - Gexp * diff(year[noNA[corr +
                                                                       0:1]])
            }
          }
          # update annual growth
          noNA = which(!is.na(dbhc))
          dG = anGrowth(dbhc[noNA], year[noNA], relat_change)
        }
      }
        return(dbhc[order(yorder)])
    }
  }
}

anGrowth = function(dbh, year, rel = FALSE) {
  if (rel) {
    dG = (dbh[-1] / dbhc[-length(dbh)]) ^ (1 / diff(year)) - 1
  } else {
    dG = diff(dbh) / diff(year)
  }
  return(dG)
}
