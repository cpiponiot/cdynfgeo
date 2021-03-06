#' Open and clean species table.
#' @param path_spp  The directory were the species table(s) are. Species tables must
#' be .Rdata files, and contain 'spp' in the name. Columns must contain: xxx
#'
#' @return A data table with the following columns: "site" (site name, as spelled in the .Rdata file), "sp" (species code), "name" (Latin name), "taxo_group" (palm, vine, tree fern, gymnosperm, or NA for other groups)
#' @export

clean_spptab = function(path_spp = getwd()) {

  list_spp = list.files(path_spp, pattern = "spp", full.names = TRUE)
  if (length(list_spp) == 0)
    stop("No species table in your directory (must contain 'spp' in the file name)")

  spptab = lapply(list_spp, open_FGEOdata)
  spptab = data.table::rbindlist(spptab, fill = TRUE)

  if (!is.null(spptab$spcode))
    spptab[!is.na(spcode)]$sp = spptab[!is.na(spcode)]$spcode

  if (!is.null(spptab$latin)) {
    nalatin = which(!is.na(spptab$genus) & !is.na(spptab$species) & is.na(spptab$latin))
    spptab[nalatin]$latin = paste(spptab[nalatin]$genus, spptab[nalatin]$species)
    spptab = subset(spptab, latin != "")
    # some genus and species are missing (but are in the latin column)
    spptab[is.na(genus)]$species = data.table::tstrsplit(spptab[is.na(genus)]$latin, " ")[[2]]
    spptab[is.na(genus)]$genus = data.table::tstrsplit(spptab[is.na(genus)]$latin, " ")[[1]]
  }
  # remove everything after species name
  spptab$species = data.table::tstrsplit(spptab$species, " ")[[1]]

  # remove undetermined species
  spptab[spptab$species %in% c("sp.", "spp.", "sp", "unkn")]$species = ""
  spptab[grep("[1-9]", species)]$species = ""
  spptab = subset(spptab,!grepl("nident", genus))

  spptab$name = paste(spptab$genus, spptab$species)
  # remove unnecessary tabs, "aff." and "cf."
  spptab$name = gsub("  ", " ", spptab$name)
  spptab$name = gsub(" $|aff.|cf.", "", spptab$name)
  # change one name with strange characters
  spptab[sp == "GONSPI"]$name = "Gonzalagunia hirsuta"

  # find species accepted names and correct typos
  taxoInfo = list()
  lname = split(unique(spptab$name), ceiling(seq_along(unique(spptab$name))/20))
  for (i in 1:length(lname)) {
    temp = try(taxize::tnrs(query = lname[[i]], source = "ncbi"), TRUE)
    if(isTRUE(class(temp)=="try-error")) { next } else {
      taxoInfo = rbind(taxoInfo, temp[, c("submittedname", "acceptedname")])
      print(paste0(i, "/", length(lname), " done"))
    }
  }
  taxoInfo$submittedname = gsub("\r", "", taxoInfo$submittedname)
  ## taxo that gave errors (are not in taxoInfo yet)
  missing = unique(spptab$name)[!unique(spptab$name) %in% taxoInfo$submittedname]
  ## 1 by 1
  for (i in 1:length(missing)) {
    temp = try(taxize::tnrs(query = missing[i], source = "ncbi"), TRUE)
    if (isTRUE(class(temp)=="try-error") | length(temp) == 0) { next } else {
      taxoInfo = rbind(taxoInfo, temp[, c("submittedname", "acceptedname")])
      print(paste0(i, "/", length(missing), " done"))
    }
  }
  taxoInfo$submittedname = gsub("\r", "", taxoInfo$submittedname)
  missing = unique(spptab$name)[!unique(spptab$name) %in% taxoInfo$submittedname]
  taxoInfo = rbind(taxoInfo, data.frame(submittedname = missing, acceptedname = missing))
  taxoInfo = unique(taxoInfo)

  spptab = merge(spptab,
                 taxoInfo,
                 by.x = "name",
                 by.y = "submittedname",
                 all.x = TRUE)
  accname = which(!is.na(spptab$acceptedname) & spptab$acceptedname != "")
  spptab[accname]$name = spptab[accname]$acceptedname

  # for cocoli and sherman, replace with bci species:
  # trusted database, no duplicates
  sp_cocoli = spptab[site == "bci"]
  sp_cocoli$site = "cocoli"
  sp_sherman = spptab[site == "bci"]
  sp_sherman$site = "sherman"
  spptab = subset(spptab,!(site %in% c("cocoli", "sherman")))
  spptab = rbind(spptab, sp_cocoli, sp_sherman)

  spptab = unique(spptab[, c("site", "sp", "name")])

  spptab$genus = data.table::tstrsplit(spptab$name, " ")[[1]]
  spptab$species = data.table::tstrsplit(spptab$name, " ")[[2]]

  # still a few duplicated names
  spptab$sp_site = paste(spptab$sp, spptab$site)
  duplsp = spptab[duplicated(sp_site), ]$sp_site
  spptab = subset(spptab,!(sp_site %in% duplsp & is.na(species)))
  spptab =  unique(spptab)[, c("site", "sp", "name", "genus")]

  ## add taxonomic groups
  data("taxoGroups")
  spptab = merge(spptab, taxoGroups[, c("genus", "taxo")], by = "genus", all.x = TRUE)

  spptab =  unique(spptab)[, c("site", "sp", "name", "taxo")]
  return(spptab)

}
