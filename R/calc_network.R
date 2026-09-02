#' DEPRECATED: Total Drainage Area
#' @description Calculates total drainage area given a dendritic
#' network and incremental areas.
#' @param x data.frame with ID, toID, and area columns.
#' @return numeric with total area.
#' @importFrom dplyr select left_join
#' @importFrom hydroloom accumulate_downstream
#' @export

calculate_total_drainage_area <- function(x) {

  warning("calculate_total_drainage_area() is deprecated. Please switch to hydroloom equivalent.")

  return(accumulate_downstream(x, "area"))

}

#' DEPRECATED: Calculate Arbolate Sum
#' @description Calculates arbolate sum given a dendritic
#' network and incremental lengths. Arbolate sum is the total length
#' of all upstream flowlines.
#' @param x data.frame with ID, toID, and length columns.
#' @return numeric with arbolate sum.
#' @export

calculate_arbolate_sum <- function(x) {

  warning("calculate_arbolate_sum is deprecated. Please switch to hydroloom equivalent.")

  return(accumulate_downstream(x, "length"))

}
