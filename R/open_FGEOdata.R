#' Open ForestGEO .Rdata census data.
#' @param file character string - name of the Rdata file to open (with its path). The name should be written
#'             as path/site.stem#census.rdata, or path/site.full#census.rdata
#'             (eg: data/bci.stem2.rdata is the 2nd census at bci, containing all stems)
#'
#' @return A data frame (or a list of data frames) with the census data.
#' @export

open_FGEOdata = function(file) {
  load(file)
  file_name_split = unlist(strsplit(tolower(file), "[[:punct:]]r|/"))
  file_name = file_name_split[length(file_name_split) - 1]

  temp <- get(file_name)

  # inconsistent uppercase use in colnames across sites -> all lower cases
  colnames(temp) = tolower(colnames(temp))

  dupl = which(duplicated(colnames(temp)))
  if (length(dupl)>0)
    colnames(temp)[dupl] = paste0(colnames(temp)[dupl], 2)

  # add site
  site <- data.table::tstrsplit(file_name, "[[:punct:]]")[[1]]
  site_col = which(colnames(temp)=="site")
  if (length(site_col)>0)
    temp = data.frame(temp)[, -site_col, with=FALSE]
  temp$site <- site

  # add census number
  census <- as.numeric(regmatches(file, gregexpr("[[:digit:]]+", file)))
  # remove previous "census" column
  census_col = which(colnames(temp)=="census")
  if (length(census_col)>0)
    temp = data.frame(temp)[, -census_col]
  temp$census <- census

  return(temp)
}
