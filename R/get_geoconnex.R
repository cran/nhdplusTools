#' discover geoconnex reference feature layers
#' @description
#' Queries the geoconnex.us reference feature server for available layers and
#' attributes.
#'
#' @return data.frame containing layers available and fields that are available to query.
#' @export
discover_geoconnex_reference <- function() {

  discover_oafeat(get("gocnx_ref_base_url", envir = nhdplusTools_env))

}

#' get geoconnex reference feature layers
#' @description
#' Queries the geoconnex reference feature server for features of interest.
#'
#' @param AOI  bbox, sf polygon or point, or a URL that will return an sf object when passed to
#' \link[sf]{read_sf}
#' @param type character the feature type chosen from \link{discover_geoconnex_reference}
#' @inheritParams query_usgs_oafeat
#' @param status boolean print status or not
#' @return sf data.frame containing requested reference features
#' @export
get_geoconnex_reference <- function(AOI,
                                    type = NULL,
                                    t_srs = NULL,
                                    buffer = 0.5,
                                    status = TRUE) {

  base <- get("gocnx_ref_base_url", envir = nhdplusTools_env)

  get_oafeat(base, AOI, type = type, t_srs = t_srs, buffer = buffer, status = status)

}

