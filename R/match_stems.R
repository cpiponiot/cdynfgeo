#' Matching unidentifed stems of one tree in ForestGEO data
#'
#' @param dbh numeric vector containing dbh values
#' @param year numeric vector (same size as dbh) containing each measure census year
#' @param acc_decr numeric: acceptable decrease between two censuses, as a proportion if relat_change is TRUE
#'             or in mm if relat_change is false
#' @param acc_incr numeric: acceptable annual increase, as a proportion (per year) if relat_change is TRUE
#'             or in mm/year if relat_change is false
#' @param relat_change logical: should the decrease and
#'        increse values be relative (proportion of total dbh) or absolute values (in mm)?
#'
#' @return numeric vector of the same length as dbh and year, with the stem ID of each measurement
#'
#' @export
#'

match_stems <- function(dbh,
                        year,
                        acc_decr = -5,
                        acc_incr = 35,
                        relat_change = FALSE) {
  # keep the order of dbh measurements
  names(dbh) = 1:length(dbh)
  
  if (anyDuplicated(year) == 0)  {
    ## only one stem
    return(rep("1", length(year)))
    
  }  else if (length(unique(year)) == 1) {
    ## only one census
    return(as.character(1:length(year)))
    
  } else {
    X = split(dbh, year)
    X = lapply(X, function(i)
      sort(i[!is.na(i)], decreasing = TRUE))
    
    ## time between two censuses
    dyear = diff(as.numeric(names(X)))
    
    # create matrix of stems dbh: each row = one stem, each column = one census
    stemDBH = matrix(NA, nrow = length(unlist(X)), ncol = length(X))
    # create a similar matrix for measurement ID
    stemN = matrix(NA, nrow = length(unlist(X)), ncol = length(X))
    
    # first column: first census stems
    stemDBH[1:length(X[[1]]), 1] = X[[1]]
    stemN[1:length(X[[1]]), 1] = names(X[[1]])
    
    # sequentially fill each column with a census dbh values
    # start with largest stem: choose the largest available stem
    # it can be associated with
    
    for (j in 2:length(X)) {
      # empty rows
      empty = apply(is.na(stemDBH), 1, all)
      
      for (i in 1:length(X[[j]])) {
        value = X[[j]][i]
        id = names(X[[j]][i])
        
        # available rows in the matrix
        available = is.na(stemDBH[, j])
        
        # select positions consistent with the previous census
        if (relat_change) {
          acc_range = c(value / (1 + acc_decr) ^ dyear[j - 1],
                        value / (1 + acc_incr) ^ dyear[j - 1])
        } else {
          acc_range = c(value - acc_decr * dyear[j - 1],
                        value - acc_incr * dyear[j - 1])
        }
        possible = (stemDBH[, j - 1] < acc_range[1] &
                      # decrease not too large
                      stemDBH[, j - 1] > acc_range[2]) # increase not too large
        
        # position finally chosen: the first one consistent with the previous census
        # and available (or an empty row)
        position = which((possible | empty) & available) [1]
        
        stemDBH[position, j] = value
        stemN[position, j] = id
      }
    }
    ### TODO check if newly recruited stems are not > (9.9 mm * acc_incr ^ dyear) and if so, send a warning
    
    Mstem = cbind(stemID = rep(1:nrow(stemN), ncol(stemN)),
                  stemN = as.numeric(c(stemN)))
    # remove NAs
    Mstem = Mstem[!is.na(Mstem[, 2]), ]
    return(as.character(Mstem[order(Mstem[, 2]), 1]))
    
  }
}
