#' Get Raindrop Trace
#' @description Uses a raindrop trace web service to trace the
#' nhdplus digital elevation model to the nearest downslope flowline.
#' @param point sfc POINT including crs as created by:
#' \code{sf::st_sfc(sf::st_point(.. ,..), crs)}
#' @param direction character \code{"up"}, \code{"down"}, or \code{"none"}.
#' Controls the portion of the split flowline that is returned along with
#' the raindrop trace line.
#' @return sf data.frame containing raindrop trace and requested
#' portion of flowline.
#' @export
get_raindrop_trace <- function(point, direction = "down") {

  point <- check_point(point)[[1]]

  url_base <- get_nldi_url(pygeo = TRUE)

  url <- paste0(url_base, "nldi-flowtrace/execution")

  allowed_direction <- c("up", "down", "none")

  if(!direction %in% allowed_direction)
    stop(paste("direction must be in",
               paste(allowed_direction, collapse = ", ")))

  return(sf_post(url, make_json_input_trace(point, direction = direction)))

}

#' Get split catchment
#' @description Uses a catchment splitting web service to retrieve
#' the portion of a catchment upstream of the point provided.
#' @param point scf POINT including crs as created by:
#' \code{sf::st_sfc(sf::st_point(.. ,..), crs)}.
#' @param upstream logical If TRUE, the entire drainage basin upstream
#' of the point provided is returned in addition to the local catchment.
#' @return sf data.frame containing the local catchment, the split portion
#' and optionally the total drainage basin.
#' @details
#' This service works within the coterminous US NHDPlusV2 domain. If the point
#' provided falls on an NHDPlusV2 flowline as retrieved from \link{get_raindrop_trace}
#' the catchment will be split across the flow line. IF the point is not
#' along the flowline a small sub catchment will typically result. As a result,
#' most users of this function will want to use \link{get_raindrop_trace} prior
#' to calls to this function.
#'
#' An attempt is made to eliminate polygon shards if they exist in the output.
#' However, there is a chance that this function will return a multipolygon
#' data.frame.
#'
#' @export
get_split_catchment <- function(point, upstream = TRUE) {

  point <- check_point(point)[[1]]

  url_base <- get_nldi_url(pygeo = TRUE)

  url <- paste0(url_base, "nldi-splitcatchment/execution")

  out <- sf_post(url, make_json_input_split(point, upstream),
                 err_mess = paste("Ensure that the point you submitted is within\n the",
                                  "coterminous US and consider trying get_raindrop_trace\ to ensure",
                                  "your point is not too close to a catchment boundary."))

  try({
    if(!is.null(out)) {
      sf::st_geometry(out) <- sf::st_sfc(lapply(sf::st_geometry(out), remove_shards),
                                         crs = sf::st_crs(out))

      if(sf::st_geometry_type(out, by_geometry = FALSE) == "GEOMETRY") {
        sf::st_geometry(out) <- sf::st_sfc(lapply(sf::st_geometry(out),
                                                  sf::st_cast, to = "MULTIPOLYGON"),
                                           crs = sf::st_crs(out))
      }
    }
  }, silent = TRUE)

  return(out)
}

remove_shards <- function(g, thresh = 0.01) {

  if(length(g) == 1) return(sf::st_polygon(g[[1]]))

  p <- sf::st_cast(sf::st_sfc(g), "POLYGON")

  a <- sf::st_area(p)

  p <- p[a > max(a) * thresh]

  if(length(p) > 1) return(sf::st_multipolygon(p))

  sf::st_polygon(p[[1]])
}

#' Get Cross Section From Point (experimental)
#' @description Uses a cross section retrieval web services to retrieve a
#' cross section given a point and specified width. Orientation is determined
#' based on direction of a the flowline found near point. This function uses
#' a 10m National Elevation Dataset request on the back end.
#' @param point sfc POINT including crs as created by:
#' \code{sf::st_sfc(sf::st_point(.. ,..), crs)}crs.
#' @param width Cross section width in meters.
#' @param num_pts numeric number of points to retrieve along the cross section.
#' @return sf data.frame containing points retrieved.
#' @export
get_xs_point <- function(point, width, num_pts) {

  point <- check_point(point)[[1]]

  url_base <- get_nldi_url(pygeo = TRUE)

  url <- paste0(url_base, "nldi-xsatpoint/execution")

  get_xs(url, make_json_input_xspt, point, width, num_pts)

}

#' Get Cross Section Endpoints (experimental)
#' @description Uses a cross section retrieval web services to retrieve a
#' cross section between two endpoints.
#' @param point1 sfc POINT including crs as created by:
#' \code{sf::st_sfc(sf::st_point(.. ,..), crs)}
#' @param point2 sfc POINT including crs.
#' @param num_pts numeric number of points to retrieve along the cross section.
#' @param res integer resolution of 3D Elevation Program data to request.
#' Must be on of: 1, 3, 5, 10, 30, 60.
#' @return sf data.frame containing points retrieved.
#' @export
get_xs_points <- function(point1, point2, num_pts, res = 1) {

  point1 <- check_point(point1)[[1]]
  point2 <- check_point(point2)[[1]]

  url_base <- get_nldi_url(pygeo = TRUE)

  url <- paste0(url_base, "nldi-xsatendpts/execution")

  check_res(res)

  get_xs(url, make_json_input_xspts, point1, point2, num_pts, res)

}

check_res <- function(res) {
  if(!res %in% c(1, 3, 5, 10, 30, 60)) {
    stop("res input must be on of 1, 3, 5, 10, 30, 60")
  }
  return(invisible(TRUE))
}

#' Get Elevation Along Path (experimental)
#' @description Uses a cross section retrieval web services to retrieve elevation
#'  along a path.
#' @param points sf data.frame containing a point column.
#' @param num_pts numeric number of points to retrieve along the cross section.
#' @param res integer resolution of 3D Elevation Program data to request.
#' Must be on of: 1, 3, 5, 10, 30, 60.
#' @param status logical
#' @return sf data.frame containing points retrieved. Names include
#' "id", "distance_m", "elevation_m", "spatial_ref", "geometry",
#' and ".group". .group tracks which input point each set of output
#' points belongs to.
#' @export
get_elev_along_path <- function(points, num_pts, res = 1, status = TRUE) {

  url_base <- get_nldi_url(pygeo = TRUE)

  url <- paste0(url_base, "nldi-xsatendpts/execution")

  check_res(res)

  get_elev(url, make_json_input_xspts, points, num_pts, res, status)

}


get_elev <- function(url, fun, points, num_pts, res, status) {

  points <- split(points, seq_len(nrow(points)))

  points <- lapply(points, check_point)

  points <- lapply(points, "[[", 1)

  data_elev <- data.frame()
  dist <- vector()

  for(i in 1:(length(points)-1)) {

    if(status)
      message(paste("Requestion segment", i, "of", (length(points)-1)))

    data <- get_xs(url, fun, points[[i]], points[[i+1]], num_pts, res)

    if(is.null(data)) {
      return(NULL)
    }

    data$.group <- i

    if(i == 1){

      if(num_pts %% 2 != 0){
        dist[[i]] <- data[[num_pts, 'distance_m']]
      } else {
        dist[[i]] <- data[[num_pts + 1, 'distance_m']]
      }

    } else {

      data[['distance_m']] <- dist[i-1] + data[['distance_m']]

      if(num_pts %% 2 != 0){
        dist[[i]] <- data[[num_pts, 'distance_m']]
      } else {
        dist[[i]] <- data[[num_pts + 1, 'distance_m']]
      }

    }

    data_elev <- rbind(data_elev, data)

  }

  data_elev
}

#' @importFrom dplyr rename
get_xs <- function(url, fun, ...) {
  sf <- sf_post(url, fun(...))

  if(is.null(sf)) {
    return(NULL)
  }

  rename(sf,
         distance_m = "distance",
         elevation_m = "elevation")
}

sf_post <- function(url, json, err_mess = "") {
  tryCatch({

    if(nhdplus_debug()) {
      message(paste(url, "\n"))
      message(json)
    }

    out <- httr::RETRY("POST", url, httr::accept_json(),
                       httr::content_type_json(),
                       body = json)

    if(out$status_code == 200) {
      sf::read_sf(rawToChar(out$content))
    } else {
      stop(rawToChar(out$content))
    }

  }, error = function(e) {
    message("Error calling processing service. \n Original error: \n", e,
            "\n", err_mess)
    NULL
  })
}


check_point <- function(p) {
  mess <- "Point must be of type sfc and have a CRS declared."

  if(inherits(p, "sf")) p <- sf::st_geometry(p)

  if(!inherits(p, "sfc")) stop(mess)

  tryCatch({
    sf::st_transform(p, 4326)
  }, error = function(e) {
    stop(paste(mess, "Original error was: \n", e))
  })
}

make_json_input_trace <- function(p, raindrop = TRUE, direction = "down") {

  jsonlite::toJSON(list(inputs = list(list(id = "lat",
                                           type = "text/plain",
                                           value = as.character(p[2])),
                                      list(id = "lon",
                                           type = "text/plain",
                                           value = as.character(p[1])),
                                      list(id = "raindroptrace",
                                           type = "text/plain",
                                           value = ifelse(raindrop,
                                                          "true", "false")),
                                      list(id = "direction",
                                           type = "text/plain",
                                           value = direction))),
                   pretty = TRUE, auto_unbox = TRUE)
}

make_json_input_split <- function(p, upstream = TRUE) {

  jsonlite::toJSON(list(inputs = list(list(id = "lat",
                                           type = "text/plain",
                                           value = as.character(p[2])),
                                      list(id = "lon",
                                           type = "text/plain",
                                           value = as.character(p[1])),
                                      list(id = "upstream",
                                           type = "text/plain",
                                           value = ifelse(upstream,
                                                          "true", "false")))),
                   pretty = TRUE, auto_unbox = TRUE)

}

make_json_input_xspt <- function(p, w, n) {

  jsonlite::toJSON(list(inputs = list(list(id = "lat",
                                           type = "text/plain",
                                           value = p[2]),
                                      list(id = "lon",
                                           type = "text/plain",
                                           value = p[1]),
                                      list(id = "width",
                                           type = "text/plain",
                                           value = w),
                                      list(id = "numpts",
                                           type = "text/plain",
                                           value = n))),
                   pretty = TRUE, auto_unbox = TRUE)
}


make_json_input_xspts <- function(p1, p2, n, r) {

  jsonlite::toJSON(list(inputs = list(list(id = "lat",
                                           type = "text/plain",
                                           value = c(p1[2], p2[2])),
                                      list(id = "lon",
                                           type = "text/plain",
                                           value = c(p1[1], p2[1])),
                                      list(id = "3dep_res",
                                           type = "text/plain",
                                           value = as.character(r)),
                                      list(id = "numpts",
                                           type = "text/plain",
                                           value = n))),
                   pretty = TRUE, auto_unbox = TRUE)
}
