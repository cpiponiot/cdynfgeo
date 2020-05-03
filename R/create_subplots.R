#' @param size subplot size, either a single numeric value or a 2 

create_subplots = function(x, y, size = 20, resample = FALSE, N = 100) {
  xgroup = as.numeric(cut(x, seq(0, round(max(x)), size)))
  ygroup = as.numeric(cut(y, seq(0, round(max(y)), size)))
  return(paste(xgroup, ygroup, sep = "_"))
  #TODO resampling
}

