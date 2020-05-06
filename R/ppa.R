#' Estimate canopy position of trees based on perfect plasticity approximation
#'
#' @param ca A numerical vector containing the crown area of indivudual trees
#'   (in ha)
#' @param area Area of the plot
#'
#' @return A vector with logical values, `TRUE` if the tree is in the canopy
#'   (according ot the PPA) and `FALSE` otherwise.
#'
#' @export
#'
ppa = function(ca, area) {
  ca_sorted = sort(ca, decreasing = TRUE)
  reorder = order(order(ca, decreasing = TRUE))
  cca = cumsum(c(0, ca_sorted[-length(ca)]))
  pcanopy = (area - cca) / ca_sorted

  if (pcanopy[1] < 0.5) {
    pcanopy[1] = 1
    warning("The biggest tree has a crown area more than double the size of the plot.")
  }

  return((pcanopy > 0.5)[reorder])
}

# ## test
#
# x = subset(fgeo_data, site=="scbi"&quadrat == "1329"&year=="2013")
# ca = x$CA
# area =0.04
