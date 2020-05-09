#' Create subplots of fixed size based on tree coordinates
#'
#' @param x x coordinate of each observation
#' @param y y coordinate of each observation
#' @param size Length of the subplot side (in m; subplots are squares).
#' @param Nmin Minimum acceptable abundance (n/ha) under which subplots are discarted. Default is 100.
#'
#' @return A vector with new subplot
#'
#' @export
#'
create_subplots = function(x, y, size = 20, Nmin = 100) {

  ## plot limits
  x0 = 20 * round(min(x / 20, na.rm = TRUE))
  x1 = 20 * round(max(x / 20, na.rm = TRUE))
  y0 = 20 * round(min(y / 20, na.rm = TRUE))
  y1 = 20 * round(max(y / 20, na.rm = TRUE))

 subplot = paste(as.numeric(cut(
    x,
    seq(x0, x1, size),
    right = FALSE,
    include.lowest = TRUE
  )),
  as.numeric(cut(
    y,
    seq(y0, y1, size),
    right = FALSE,
    include.lowest = TRUE
  )),
  sep = "_")

  ## trees outside
  subplot[grepl("NA", subplot)] = NA
  # quadrats with fewer then 6 stems in all years: not real quadrats
  Nsubplot = data.frame(table(subplot))
  not_subplot = Nsubplot$subplot[Nsubplot$Freq/((size/100)**2) < Nmin]
  subplot[subplot %in% not_subplot] = NA

  return(subplot)
}

