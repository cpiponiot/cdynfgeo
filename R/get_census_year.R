#' Extract census year from a set of exact measure dates
#' @param date vector of exact measurement dates, containing the year,
#'             month and day (not always in this order) and sometimes
#'             the hour of all measurements during one census
#'
get_census_year = function(date) {
  
  # remove hours
  split_date = tstrsplit(date, " ")[[1]]
  
  # separate year, month and day
  split_date = tstrsplit(split_date, "-|/")
  
  # the one in the middle cannot be the year
  split_date = list(split_date[[1]], split_date[[3]])
  
  split_date = lapply(split_date, as.numeric)
  
  ## remove empty dates (zeros)
  zeros = which(split_date[[1]]==0 & split_date[[2]]==0)
  if (length(zeros) > 0)
    split_date = lapply(split_date, function(x) x[-zeros])
  
  # because they are not always in the same order (depending on the site),
  # and the year is in different formats (2 or 4 digits)
  # we chose the value > 31 when there is one,
  # or the value with the smallest range across tree measurements
  
  if (any(split_date[[1]] > 31, na.rm = TRUE)) {
    year = split_date[[1]]
  } else if (any(split_date[[2]] > 31, na.rm = TRUE)) {
    year = split_date[[2]]
  } else {
    range1 = diff(range(split_date[[1]], na.rm = TRUE))
    range2 = diff(range(split_date[[2]], na.rm = TRUE))
    
    if ((range1 > 2 & range2 > 2) | range1 == range2)
      stop("There is a problem with the date format.")
    
    year = split_date[[which.min(c(range1, range2))]]
  }
  
  ## if census year has less than 4 digits: add 19- or 20-
  if (all(year < 100, na.rm = TRUE))   {
    add_1900 = which(!is.na(year) & year >= 50)
    add_2000 = which(!is.na(year) & year < 50)
    year[add_1900] = 1900 + year[add_1900]
    year[add_2000] = 2000 + year[add_2000]
  }
  
  census_year = round(mean(year, na.rm = TRUE))
  
  
  return(census_year)
  
}