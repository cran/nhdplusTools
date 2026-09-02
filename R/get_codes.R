#' @title Get Streamorder
#' @description Applies a topological sort and calculates strahler stream order.
#' Algorithm: If more than one upstream flowpath has an order equal to the
#' maximum upstream order then the downstream flowpath is assigned the maximum
#' upstream order plus one. Otherwise it is assigned the max upstream order.
#' @param x data.frame with dendritic ID and toID columns.
#' @param status logical show progress update messages?
#' @return numeric stream order in same order as input
#' @importFrom dplyr left_join select
#' @importFrom hydroloom add_streamorder
#' @export
get_streamorder <- function(x, status = TRUE) {
  check_names(x, "get_streamorder")

  add_streamorder(x, status)$stream_order

}

#' @title Get Streamlevel
#' @description Applies a topological sort and calculates stream level.
#' Algorithm: Terminal level paths are assigned level 1 (see note 1).
#' Paths that terminate at a level 1 are assigned level 2. This pattern is
#' repeated until no paths remain.
#'
#' If a TRUE/FALSE coastal attribute is included, coastal terminal paths
#' begin at 1 and internal terminal paths begin at 4 as is implemented by
#' the NHD stream leveling rules.
#'
#' @param x data.frame with levelpathi, dnlevelpat, and optionally a
#' coastal flag. If no coastal flag is included, all terminal paths are
#' assumed to be coastal.
#'
#' @return numeric stream order in same order as input
#' @importFrom hydroloom add_streamlevel
#' @export
get_streamlevel <- function(x) {

  check_names(x, "get_streamlevel")

  coastal <- NULL
  coastal <- if("coastal" %in% names(x)) "coastal"

  # TODO: drop suppressWarnings once hydroloom's hy() stops probing $id
  # for streamlevel-only inputs (warns "Unknown column: id").
  suppressWarnings(add_streamlevel(x, coastal))$stream_level

}

#' @title Get Pfafstetter Codes (DEPRECATED)
#' @description Determines Pfafstetter codes for a dendritic network with
#' total drainage area, levelpath, and topo_sort attributes.
#' @param x sf data.frame with ID, toID, totda, outletID, topo_sort,
#' and levelpath attributes.
#' @param max_level integer number of pfaf levels to attempt to calculate.
#' If the network doesn't have resolution to support the desired level,
#' unexpected behavior may occur.
#' @param status boolean print status or not
#' @importFrom hydroloom add_pfafstetter
#' @return data.frame with ID and pfaf columns.
#' @export
get_pfaf <- function(x, max_level = 2, status = FALSE) {

  warning("get_pfaf is deprecated, please use hydroloom")

  x <- st_drop_geometry(x)

  check_names(x, "get_pfaf")

  x <- select(x, all_of(get("get_pfaf_attributes",
                            envir = nhdplusTools_env)))

  x <- add_pfafstetter(x, max_level, status)

  x <- select(x, -all_of(c("toID", "totda", "outletID", "topo_sort", "levelpath")))

  x <- filter(x, !is.na(.data$pf_level_1))

  return(x)
}
