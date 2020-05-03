#' Visualize corrections made by the consolidate_data() function
#'
#' @param df data.table: The output of the function consolidate_data()
#' @param dbhmin numeric: the minimum dbh of trees to plot 
#' @param plotsPerPage numeric: number of plots (stemids) per page
#'
#' @return A list of ggplots
#'
#' @export

visualize_corrections = function(df, dbhmin = 500, plotsPerPage = 25) {
  
  library(ggplot2)
  
  big_corr = unique(df[dbhcc != dbhtc & dbhc > dbhmin]$stemid)
  dfcorr = df[stemid %in% big_corr]
  dfcorr = melt(dfcorr, id.vars = c("stemid", "year", "name"), 
                measure.vars = c("dbh", "dbhtc", "dbhcc"), 
                value.name = "dbh")
  levels(dfcorr$variable) = c("original", "taper", "final")
  final = dfcorr[variable == "final"]
  dfcorr = dfcorr[!duplicated(dfcorr[,-"variable"])]
  
  library(ggplot2)
  library(ggpubr)
  
  Np = ceiling(length(big_corr)/plotsPerPage) # number of pages
  groups = rep(1:Np, each = plotsPerPage)[1:length(big_corr)]
  pages = split(big_corr, groups)
  
  plt_page = function(ids) {
    ggplot(dfcorr[stemid %in% ids], aes(x=year, y = dbh, color = variable)) + 
      geom_point() + 
      geom_line(data = dfcorr[stemid %in% ids & variable == "original"]) +
      geom_line(data = final[stemid %in% ids]) +
      # scale_color_manual(values = c("red","darkblue", "forestgreen")) +
      facet_wrap( ~ paste(stemid, name, sep = "\n"), scales = "free")
  }
  
  plots = lapply(pages, plt_page)
  return(plots)
}