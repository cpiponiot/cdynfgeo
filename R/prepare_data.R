#' Import and homogenize ForestGEO-like data.
#'
#' Main routine to format, detect major obvious errors, and gap-fill those
#' errors in ForestGEO-like data.
#'
#' @param path String giving a path to a parent directory containing species
#'   and census datasets.
#' @param stem `TRUE`(default) or `FALSE` reflect that your censuses are
#'   ForestGEO `stem` or `tree` (aka `full`) tables, respectively.
#' @param site String giving the name of your site -- one of one of
#'   `site_info$site` (e.g., 'barro colorado island'.
#' @param dbh_units set the unit ("mm" or "cm") of DBH values, by default
#'   dbh_units=="mm".
#' @param exclude_interval `NULL` by default. If needed a vector (e.g. c(1,2))
#'   indicating which census interval(s) must be discarded from computation due,
#'   for instance, to a change in measurement protocol.
#' @param keep_column character vector containing the name of the columns to keep
#'   in the data. Default is: c("treeid", "stemid", "tag", "stemtag",
#'                              "sp", "quadrat", "gx", "gy", "dbh",
#'                              "hom", "exactdate", "codes", "dfstatus",
#'                              "agb", "brokenladder", "fire.code")
#'
#' @return A data.table (data.frame) with all relevant variables.
#'
#' @export
#'
prepare_data <- function(path,
                         site,
                         stem = TRUE,
                         dbh_units = "mm",
                         exclude_interval = NULL,
                         keep_columns = c("treeid", "stemid", "tag", "stemtag",
                                          "sp", "quadrat", "gx", "gy", "dbh",
                                          "hom", "exactdate", "codes", "dfstatus",
                                          "agb", "brokenladder", "fire.code") # for Yosemite
)
{

  print(paste("Preparing", site, "data."))
  data("site_info")

  # check that the site is the list of ForestGEO sites
  site <- tolower(site)
  INDEX <- match(tolower(site), site_info$site)
  if (is.na(INDEX)) {
    stop("Site name should be one of the following: \n",
         paste(levels(factor(site_info$site)), collapse = " - "),
         call. = FALSE)
  }

  # open census data
  path_folder <- getwd()
  files <- list.files(path, pattern = site)

  ifelse(stem,
         files <- files[grep("stem", files)],
         files <- files[grep("full", files)])

  # Check that there are files with site.stem or site.full names
  if (length(files) == 0) {
    ifelse(stem,
           stop("There are no files with site.stem names."),
           stop("There are no files with site.full names."))
  }

  # Check that all selected files are in the Rdata format
  if (!all(grepl(".rda", tolower(files))))
    stop("Census data should be in the Rdata format (.rda or .rdata extension).")

  # search for census number in file name
  census_number <-
    as.numeric(regmatches(files, gregexpr("[[:digit:]]+", files)))

  # Remove censuses to exclude (if needed)
  if (!is.null(exclude_interval) &
      sum(exclude_interval %in% census_number) > 0) {
    files <- files[-which(census_number == exclude_interval)]
    print(paste("Census", exclude_interval, "excluded."))
  }

  # open data
  census_list <- lapply(paste(path, files, sep = "/"), open_FGEOdata)
  print(paste("Opened", site, "census data."))

  DT <- data.table::rbindlist(census_list, fill = TRUE)

  # numeric columns
  DT$dbh = as.numeric(as.character(DT$dbh))
  if (is.null(DT$gx) & !is.null(DT$px)) {
    DT$gx = DT$px
    DT$gy = DT$py
  }
  DT$gx = as.numeric(as.character(DT$gx))
  DT$gy = as.numeric(as.character(DT$gy))

  # if dbh is in cm: change to mm
  if (dbh_units == "cm" |
      (!is.na(site_info$dbh_units[INDEX]) &
       site_info$dbh_units[INDEX] == "cm")) {
    DT$dbh = 10 * DT$dbh
    print("Converted dbh from cm to mm.")
  }

  # get mean year of census from exactdate
  if (is.null(DT$exactdate) & length(grep("date", colnames(DT))) == 1)
    DT$exactdate = data.frame(DT)[, grep("date", colnames(DT))]
  DT$exactdate = as.character(DT$exactdate) ## convert to always have characters

  # calculate census year (per quadrat if info is available)
  if (!is.null(DT$quadrat)) {
    DT$quadrat = as.character(DT$quadrat)
    cq_year = tapply(DT$exactdate, paste(DT$census, DT$quadrat), get_census_year)
    cq_year = data.frame(census_quadrat = names(cq_year), year = cq_year)
    cq_year$census = as.numeric(data.table::tstrsplit(cq_year$census_quadrat, " ")[[1]])
    cq_year$quadrat = data.table::tstrsplit(cq_year$census_quadrat, " ")[[2]]
    ## complete for quadrat with missing dates:
    census_year = tapply(DT$exactdate, DT$census, get_census_year)
    census_year = data.frame(census = as.numeric(names(census_year)), mean_year = census_year)
    cq_year = merge(cq_year, census_year, by = "census")
    cq_year$year[is.na(cq_year$year)] = cq_year$mean_year[is.na(cq_year$year)]
    DT = merge(DT, cq_year[, c("census", "quadrat", "year")], by = c("census", "quadrat"))
  } else  {
    census_year = tapply(DT$exactdate, DT$census, get_census_year)
    census_year = data.frame(census = as.numeric(names(census_year)), year = census_year)
    DT = merge(DT, census_year, by = "census")
  }

  # fill empty DFstatus
  if (is.null(DT$dfstatus) & length(grep("status", colnames(DT))) == 1)
    DT$dfstatus = data.frame(DT)[, grep("status", colnames(DT))]
  if (length(grep("status", colnames(DT)))==0)
    DT$dfstatus = c("alive","dead")[is.na(DT$dbh)+1]

  # transform tag number into unique tree ID and stem ID
  if (length(DT$tag) >  0 & length(DT$treeid) == 0)
    DT$treeid = as.numeric(as.factor(DT$tag))
  if (length(DT$treeid) > 0)
    DT$treeid = paste(site, DT$treeid, sep = "_")
  if (length(DT$stemtag) >  0 & length(DT$stemid) == 0)
    DT$stemid = as.numeric(as.factor(DT$stemtag))
  if (length(DT$stemid) > 0)
    DT$stemid = paste(DT$treeid, DT$stemid, sep = "_")

  # pom and hom are the same measure: replace missing hom by pom
  if (all(c("hom", "pom") %in% colnames(DT))) {
    DT$hom = as.numeric(DT$hom)
    DT[is.na(hom)]$hom = as.numeric(DT[is.na(hom)]$pom)
  } else if ("pom" %in% colnames(DT)) {
    DT$hom = as.numeric(DT$pom)
  } else {
    DT$hom = as.numeric(DT$hom)
  }

  # spcode instead of sp
  if (sum(colnames(DT)=='sp')==0 & length(DT$spcode) > 0)
    DT$sp = DT$spcode
  if (sum(colnames(DT)=='sp')==0 & length(DT$mnemonic) > 0)
    DT$sp = DT$mnemonic

  # change all "" with NA in stem tags
  DT$stemtag[DT$stemtag == ""] = NA

  keep_columns = c("site", "census", "year", keep_columns)
  DT = DT[, which(colnames(DT) %in% keep_columns), with = FALSE]

  return(DT)
}

#' Extract census year from a set of exact measure dates
#' @param date vector of exact measurement dates, containing the year,
#'             month and day (not always in this order) and sometimes
#'             the hour of all measurements during one census
#' @return A numeric vector with census years
#' @export
#'
