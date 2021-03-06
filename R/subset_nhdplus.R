#' @title Subset NHDPlus
#' @description Saves a subset of the National Seamless database or other
#' nhdplusTools compatible data based on a specified collection of COMIDs.
#' This function uses \code{\link{get_nhdplus}} for the "download" data
#' source but returns data consistent with local data subsets in a subset
#' file.
#' @param comids integer vector of COMIDs to include.
#' @param output_file character path to save the output to defaults
#' to the directory of the nhdplus_data.
#' @param nhdplus_data character path to the .gpkg or .gdb containing
#' the national seamless database, a subset of NHDPlusHR,
#' or "download" to use a web service to download NHDPlusV2.1 data.
#' Not required if \code{\link{nhdplus_path}} has been set or the default
#' has been adopted. See details for more.
#' @param bbox object of class "bbox" as returned by sf::st_bbox in Latitude/Longitude.
#' If no CRS is present, will be assumed to be in WGS84 Latitude Longitude.
#' @param simplified boolean if TRUE (the default) the CatchmentSP layer
#' will be included. Not relevant to the "download" option or NHDPlusHR data.
#' @param overwrite boolean should the output file be overwritten
#' @param return_data boolean if FALSE path to output file is returned silently otherwise
#' data is returned in a list.
#' @param status boolean should the function print status messages
#' @param flowline_only boolean WARNING: experimental
#' if TRUE only the flowline network and attributes will be returned
#' @param streamorder integer only streams of order greater than or equal will be downloaded.
#' Not implemented for local data.
#' @param out_prj character override the default output CRS of NAD83 lat/lon (EPSG:4269)
#' @details
#'
#' This function relies on the National Seamless Geodatabase or Geopackage.
#' It can be downloaded
#' \href{https://www.epa.gov/waterdata/nhdplus-national-data}{here.}
#'
#' The "download" option of this function should be considered preliminary
#' and subject to revision. It does not include as many layers and may not
#' be available permenently.
#'
#' @return character path to the saved subset geopackage
#' @export
#' @examples
#' \donttest{
#'
#' source(system.file("extdata/sample_data.R", package = "nhdplusTools"))
#'
#' nhdplus_path(sample_data)
#'
#' staged_nhdplus <- stage_national_data(output_path = tempdir())
#'
#' sample_flines <- readRDS(staged_nhdplus$flowline)
#'
#' geom_col <- attr(sample_flines, "sf_column")
#'
#' plot(sample_flines[[geom_col]],
#'      lwd = 3)
#'
#' start_point <- sf::st_sfc(sf::st_point(c(-89.362239, 43.090266)),
#'                           crs = 4326)
#'
#' plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
#'
#' start_comid <- discover_nhdplus_id(start_point)
#'
#' comids <- get_UT(sample_flines, start_comid)
#'
#' plot(dplyr::filter(sample_flines, COMID %in% comids)[[geom_col]],
#'      add=TRUE, col = "red", lwd = 2)
#'
#' output_file <- tempfile(fileext = ".gpkg")
#'
#' subset_nhdplus(comids = comids,
#'                output_file = output_file,
#'                nhdplus_data = sample_data,
#'                overwrite = TRUE,
#'                status = TRUE)
#'
#' sf::st_layers(output_file)
#'
#' catchment <- sf::read_sf(output_file, "CatchmentSP")
#'
#' plot(catchment[[attr(catchment, "sf_column")]], add = TRUE)
#'
#' waterbody <- sf::read_sf(output_file, "NHDWaterbody")
#'
#' plot(waterbody[[attr(waterbody, "sf_column")]],
#'      col = rgb(0, 0, 1, alpha = 0.5), add = TRUE)
#'
#' # Cleanup temp
#' sapply(staged_nhdplus, unlink)
#' unlink(output_file)
#'
#' # Download Option:
#' subset_nhdplus(comids = comids,
#'                output_file = output_file,
#'                nhdplus_data = "download",
#'                overwrite = TRUE,
#'                status = TRUE, flowline_only = FALSE)
#'
#' sf::st_layers(output_file)
#'
#' # NHDPlusHR
#' source(system.file("extdata/nhdplushr_data.R", package = "nhdplusTools"))
#'
#' up_ids <- get_UT(hr_data$NHDFlowline, 15000500028335)
#'
#' sub_gpkg <- file.path(work_dir, "sub.gpkg")
#' sub_nhdhr <- subset_nhdplus(up_ids, output_file = sub_gpkg,
#'                             nhdplus_data = hr_gpkg, overwrite = TRUE)
#'
#' sf::st_layers(sub_gpkg)
#' names(sub_nhdhr)
#'
#' plot(sf::st_geometry(hr_data$NHDFlowline), lwd = 0.5)
#' plot(sf::st_geometry(sub_nhdhr$NHDFlowline), lwd = 0.6, col = "red", add = TRUE)
#'
#' unlink(output_file)
#' unlink(sub_gpkg)
#' }
#'

subset_nhdplus <- function(comids = NULL, output_file = NULL, nhdplus_data = NULL, bbox = NULL,
                           simplified = TRUE, overwrite = FALSE, return_data = TRUE, status = TRUE,
                           flowline_only = NULL, streamorder = NULL, out_prj = 4269) {

  if(is.null(flowline_only)) {
    if(!is.null(nhdplus_data) && nhdplus_data == "download") {
      flowline_only <- TRUE
    } else {
      flowline_only <- FALSE
    }
  }

  if(!is.null(comids) && length(comids) == 0) stop("comids must be NULL or non-empty")

  if (status) message("All intersections performed in latitude/longitude.")

  if(any(bbox > 180 | bbox < -180)) stop("invalid bbox entry")

  if(!is.null(bbox) && nhdplus_data == "download") {
    x_range <- bbox[3] - bbox[1]
    y_range <- bbox[4] - bbox[2]

    if((x_range * y_range) > 10) {
      warning("Large bounding box submitted for download. Performance may be slow or unstable.")
    }
  }

  if(!is.null(output_file)) {
    if (!grepl("*.gpkg$", output_file)) {
      stop("output_file must end in '.gpkg'")
    }

    if (file.exists(output_file) & !overwrite) {
      stop("output_file exists and overwrite is false.")
    } else if (file.exists(output_file) & overwrite) {
      unlink(output_file)
    }
  }

  if (is.null(nhdplus_data)) {
    nhdplus_data <- nhdplus_path()
  }

  check_nhd_data(nhdplus_data)

  if(is.null(bbox)) {
    if(is.null(comids)) stop("must provide comids or bounding box")

    out_list <- c(get_flowline_subset(nhdplus_data, comids,
                                      output_file,
                                      status, out_prj))

    if(!flowline_only) {
      tryCatch({
        out_list <- c(out_list, get_catchment_subset(nhdplus_data, comids,
                                                     output_file, simplified,
                                                     status, out_prj))

        catch_layer <- get_catchment_layer_name(simplified, nhdplus_data)

        envelope <- sf::st_transform(sf::st_as_sfc(sf::st_bbox(out_list[[catch_layer]])),
                                     4326)

        intersection_names <- c("NHDArea", "NHDWaterbody", "NHDFlowline_NonNetwork")
      }, error = function(e) {
        warning(paste("error getting catchment from nhdplus_data\n", e))
      })

      if(!exists("intersection_names")) intersection_names <- c()

    } else {
      intersection_names <- c()
    }

  } else {
    out_list <- list()

    if(!is.null(comids)) warning("using bounding box rather than submitted comids")

    if(!is.null(attr(bbox, "crs"))) {
      envelope <- sf::st_transform(sf::st_as_sfc(bbox),
                                   4326)
    } else {
      if((length(bbox) != 4 | !is.numeric(bbox)) |
         (!(all(bbox >= -180) & all(bbox <= 180)))) stop("invalid bbox entry")
      names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
      bbox <- sf::st_bbox(bbox, crs = sf::st_crs(4326))
      envelope <- sf::st_as_sfc(bbox)
    }

    intersection_names <- c(get_catchment_layer_name(simplified, nhdplus_data),
                            get_flowline_layer_name(nhdplus_data),
                            "NHDArea", "NHDWaterbody", "NHDFlowline_NonNetwork")

    if(flowline_only) intersection_names <- get_flowline_layer_name(nhdplus_data)
  }

  if (nhdplus_data == "download") {

    for (layer_name in intersection_names) {
      layer <- sf::st_transform(envelope, 4326) %>%
        get_nhdplus_bybox(layer = tolower(layer_name), streamorder = streamorder)

      if(!is.null(nrow(layer)) && nrow(layer) > 0) {
        layer <- check_valid(layer, out_prj)

        if(return_data) {
          out_list[layer_name] <- list(layer)
        }

        if(!is.null(output_file)) {
          sf::write_sf(clean_bbox(layer), output_file, layer_name)
        }
      }
    }

  } else {
    if(!flowline_only) {
      if("Gage" %in% st_layers(nhdplus_data)$name) {
        intersection_names <- c(intersection_names, "Gage", "Sink")
      } else {
        intersection_names <- c(intersection_names, "NHDPlusSink")
        intersection_names <- intersection_names[which(intersection_names %in% st_layers(nhdplus_data)$name)]
      }
    }

    out_list <- c(out_list,
                  stats::setNames(lapply(intersection_names, intersection_write,
                                         data_path = nhdplus_data,
                                         envelope = envelope,
                                         output_file = output_file,
                                         status = status,
                                         out_prj = out_prj), intersection_names))
  }

  if(return_data) return(out_list)

  return(output_file)
}

intersection_write <- function(layer_name, data_path, envelope,
                               output_file, status, out_prj) {
  out_list <- list()

  if (status) message(paste("Reading", layer_name))
  layer <- sf::st_zm(sf::read_sf(data_path, layer_name))

  intersection_test <- c()

  try(intersection_test <- suppressMessages(sf::st_intersects(
    sf::st_transform(layer, 4326), envelope)), silent = TRUE)

  found <- lengths(intersection_test)

  if(length(found) > 0) {
    out <- dplyr::filter(layer, found > 0)
  } else {
    out <- data.frame()
  }

  if (nrow(out) > 0) {

    out <- check_valid(out, out_prj)

    if (status) message(paste("Writing", layer_name))
    if (is.null(output_file)) {
      return(out)
    } else {
      sf::write_sf(clean_bbox(out), output_file, layer_name)
      return(layer_name)
    }
  } else {
    if (status) message(paste("No features to write in", layer_name))
  }
}

#' @title Stage NHDPlus National Data (deprecated)
#' @description Breaks down the national geo database into a collection
#' of quick to access R binary files.
#' @param include character vector containing one or more of:
#' "attributes", "flowline", "catchment".
#' @param output_path character path to save the output to defaults
#' to the directory of the nhdplus_data.
#' @param nhdplus_data character path to the .gpkg or .gdb
#' containing the national seamless dataset. Not required if
#' \code{\link{nhdplus_path}} has been set.
#' @param simplified boolean if TRUE (the default) the CatchmentSP layer
#' will be included.
#' @details "attributes" will save `NHDFlowline_Network` attributes
#' as a separate data.frame without the geometry. The others will save
#' the `NHDFlowline_Network` and `Catchment` or `CatchmentSP`
#' (per the `simplified` parameter) as sf data.frames with
#' superfluous Z information dropped.
#'
#' The returned list of paths is also added to the nhdplusTools_env
#' as "national_data".
#'
#' @return list containing paths to the .rds files.
#' @export
#' @examples
#'
#' source(system.file("extdata/sample_data.R", package = "nhdplusTools"))
#'
#' stage_national_data(nhdplus_data = sample_data, output_path = tempdir())
#'
stage_national_data <- function(include = c("attribute",
                                            "flowline",
                                            "catchment"),
                                output_path = NULL,
                                nhdplus_data = NULL,
                                simplified = TRUE) {

  if (is.null(output_path)) {
    output_path <- dirname(nhdplus_path())
    warning(paste("No output path provided, using:", output_path))
  }

  if (is.null(nhdplus_data)) {
    nhdplus_data <- nhdplus_path()

    if (nhdplus_data == get("default_nhdplus_path",
                            envir = nhdplusTools_env) &
        !file.exists(nhdplus_data)) {
      stop(paste("Didn't find NHDPlus national data in default location:",
                 nhdplus_data))
    } else if (!file.exists(nhdplus_data)) {
      stop(paste("Didn't find NHDPlus national data",
                 "in user specified location:",
                 nhdplus_data))
    }
  }

  allow_include <- c("attribute", "flowline", "catchment")

  if (!all(include %in% allow_include)) {
    stop(paste0("Got invalid include entries. Expect one or more of: ",
                paste(allow_include, collapse = ", "), "."))
  }

  outlist <- list()

  if (any(c("attribute", "flowline") %in% include)) {

    out_path_attributes <- file.path(output_path,
                                     "nhdplus_flowline_attributes.rds")
    out_path_flines <- file.path(output_path, "nhdplus_flowline.rds")

    if (!(file.exists(out_path_flines) | file.exists(out_path_attributes))) {
      fline <- sf::st_zm(sf::read_sf(nhdplus_data,
                                     get_flowline_layer_name(nhdplus_data)))
    }

    if ("attribute" %in% include) {
      if (file.exists(out_path_attributes)) {
        warning("attributes file exists")
      } else {
        saveRDS(sf::st_set_geometry(fline, NULL), out_path_attributes)
      }
      outlist["attributes"] <- out_path_attributes
    }

    if ("flowline" %in% include) {
      if (file.exists(out_path_flines)) {
        warning("flowline file exists")
      } else {
        saveRDS(fline, out_path_flines)
      }
      outlist["flowline"] <- out_path_flines
    }
  }

  if (exists("fline")) rm(fline)

  if ("catchment" %in% include) {
    out_path_catchments <- file.path(output_path, "nhdplus_catchment.rds")
    if (file.exists(out_path_catchments)) {
      warning("catchment already exists.")
    } else {

      layer_name <- get_catchment_layer_name(simplified, nhdplus_data)

      saveRDS(sf::st_zm(sf::read_sf(nhdplus_data, layer_name)),
              out_path_catchments)
    }
    outlist["catchment"] <- out_path_catchments
  }
  assign("national_data", outlist, envir = nhdplusTools_env)

  return(outlist)
}

#' @title Try to find staged NHDPlus data
#' @noRd
check_nhd_data <- function(nhdplus_data) {

  if(nhdplus_data == "download") return(invisible(TRUE))

  if(file.exists(nhdplus_data)) {
    return(invisible(TRUE))
  } else {
    stop("couldn't find nhdplus data")
  }
}

#' @title Get subset of flowline data later.
#' @noRd
get_flowline_subset <- function(nhdplus_data, comids, output_file,
                                status, out_prj) {

  layer_name <- get_flowline_layer_name(nhdplus_data)

  if (status) message(paste("Reading", layer_name))

  if (nhdplus_data == "download") {

    if (length(comids) > 3000) {
      warning("Download functionality not tested for this many comids")
    }

    fline <- get_nhdplus_byid(comids, tolower(layer_name))

  } else {

    if(!layer_name %in% st_layers(nhdplus_data)$name) {
      layer_name <- "NHDFlowline"
    }

    fline <- get_nhd_data(nhdplus_data,layer_name, comids, "COMID", status)

  }

  fline <- check_valid(fline, out_prj)

  if (status) message(paste("Writing", layer_name))

  if(!is.null(output_file)) {
    sf::write_sf(clean_bbox(fline), output_file, layer_name)
  }
  out <- list()
  out[layer_name] <- list(fline)

  return(out)
}

get_nhd_data <- function(nhdplus_data, layer_name, comids, id, status) {

  sets <- lapply(1:ceiling(length(comids) / 1000), function(x) {
    start <- 1000 * (x - 1) + 1

    end <- 1000 * x

    end <- ifelse(end > length(comids), length(comids), end)

    comids[start:end]
  })

  assign("cur_count", 0, envir = nhdplusTools_env)

  out <- lapply(sets, function(x, total) {
    if(status) {

      cur_count <-
        get("cur_count", envir = nhdplusTools_env) + length(x)

      assign("cur_count", cur_count, envir = nhdplusTools_env)

      message(paste(cur_count, "comids of", total))

    }

    align_nhdplus_names(
      sf::read_sf(nhdplus_data, layer_name,
                  query = get_query(nhdplus_data, layer_name,
                                    id, x)))
  }, total = sum(lengths(sets)))

  do.call(rbind, out)

}

get_query <- function(nhdplus_data, layer_name, id, comids) {
  layer_atts <- sf::read_sf(nhdplus_data, layer_name,
                            query = paste0("SELECT * FROM ",
                                           layer_name, " LIMIT 1"))

  update_atts <- align_nhdplus_names(layer_atts)

  id_att <- names(layer_atts)[which(names(update_atts) == id)]

  query <- paste0("SELECT * FROM ", layer_name,
                  " WHERE ", id_att, " IN (",
                  paste(comids, collapse = ", "), ")")
}

#' @title Get subset of catchment data layer.
#' @noRd
get_catchment_subset <- function(nhdplus_data, comids, output_file,
                                 simplified, status, out_prj) {

  layer_name <- get_catchment_layer_name(simplified, nhdplus_data)

  if (status) message(paste("Reading", layer_name))

  if (nhdplus_data == "download") {

    catchment <- get_nhdplus_byid(comids, tolower(layer_name))

  } else {

    catchment <- get_nhd_data(nhdplus_data,layer_name, comids, "FEATUREID", status)

  }

  catchment <- check_valid(catchment, out_prj)

  if (status) message(paste("Writing", layer_name))

  if(!is.null(output_file)) {
    sf::write_sf(clean_bbox(catchment), output_file, layer_name)
  }
  out <- list()
  out[layer_name] <- list(catchment)
  return(out)
}

clean_bbox <- function(x) {
  if("bbox" %in% names(x) && class(x$bbox[1]) == "list") {
    x$bbox <- sapply(x$bbox, paste, collapse = ",")
  }

  return(x)
}

fix_g_type <- function(g, type = "POLYGON", orig_type = "MULTIPOLYGON") {

  tryCatch({
    sf::st_cast(sf::st_sfc(g[grepl(type, sapply(g, sf::st_geometry_type))]),
                orig_type)
  }, error = function(e) {
    sf::st_sfc(g)
  })

}

check_valid <- function(x, out_prj = sf::st_crs(x)) {

  if(is.null(x)){return(NULL)}

  x <- sf::st_zm(x)

  if (!all(sf::st_is_valid(x))) {

    message("Found invalid geometry, attempting to fix.")

    orig_type <- unique(as.character(sf::st_geometry_type(x)))

    orig_type <- orig_type[grepl("POLY|LINE", orig_type)]

    try({
      x <- sf::st_make_valid(x)

      if(!all(sf::st_geometry_type(x) == orig_type)) {
        if(any(grepl("^GEOMETRY", sf::st_geometry_type(x)))) {

          sf::st_geometry(x) <-
            sf::st_sfc(do.call(c, lapply(sf::st_geometry(x), fix_g_type,
                                         type = gsub("^MULTI", "", orig_type))),
                       crs = sf::st_crs(x))

          x <- sf::st_cast(x, orig_type)

        }
      }
    })
  }

  if (any(grepl("POLYGON", class(sf::st_geometry(x))))) {
    suppressMessages(suppressWarnings(x <- sf::st_buffer(x, 0)))
  }

  if (sf::st_crs(x) != sf::st_crs(out_prj)) {
    x <- sf::st_transform(x, out_prj)
  }

  types <- as.character(sf::st_geometry_type(x, by_geometry = FALSE))

  if(grepl("^GEOME", types)) {
    unq <- unique(as.character(
      sf::st_geometry_type(x, by_geometry = TRUE)))

    cast_to <- unq[which.max(tabulate(match(types, unq)))]

    if(any(grepl("^MULTI", unq)) & !grepl("^MULTI", cast_to)) {
      cast_to <- paste0("MULTI", cast_to)
    }

    tryCatch(x <- suppressWarnings(sf::st_cast(x, cast_to)),
             error = function(e) {
               warning(paste0("\n\n Failed to unify output geometry type. \n\n",
                             e,
                             "\n Dropping non-", cast_to, " geometries. \n"))
             })

    r <- nrow(x)

    x <- x[sf::st_geometry_type(x, by_geometry = TRUE) == cast_to, ]

    if(r != nrow(x)) {
      x <- sf::st_cast(x, cast_to)
    }
  }

  suppressWarnings(x <- sf::st_simplify(x, dTolerance = 0))

  return(x)
}

get_catchment_layer_name <- function(simplified, nhdplus_data) {
  if(is.null(nhdplus_data) || nhdplus_data == "download") { # Only simplified via download
    layer_name <- "CatchmentSP"
  } else {
    if(simplified) { # Can get simplified from local data
      layer_name <- "CatchmentSP"
    } else { # Has to be full catchment
      layer_name <- "Catchment"
    }
    if(!layer_name %in% sf::st_layers(nhdplus_data)$name) # unless it's high res.
      layer_name <- "NHDPlusCatchment"
  }
  return(layer_name)
}

get_flowline_layer_name <- function(nhdplus_data) {
  layer_name <- "NHDFlowline_Network"
  if(nhdplus_data != "download" &&
     !is.null(nhdplus_data) &&
     !layer_name %in% sf::st_layers(nhdplus_data)$name) # unless it's high res.
    layer_name <- "NHDFlowline"
  layer_name
}

#' Subset by Raster Processing Unit.
#' @description Given flowlines and an rpu_code, performs a network-safe subset such
#' that the result can be used in downstream processing. Has been tested to work
#' against the entire NHDPlusV2 domain and satisfies a number of edge cases.
#' @param fline sf data.frame NHD Flowlines with COMID, Pathlength, LENGTHKM, and Hydroseq.
#' LevelPathI, RPUID, ToNode, FromNode, and ArbolateSu.
#' @param rpu character e.g. "01a"
#' @param run_make_standalone boolean should the run_make_standalone function be run on result?
#' @export
#' @return data.frame containing subset network
#' @importFrom dplyr filter arrange summarize
#' @importFrom sf st_sf st_drop_geometry
#' @examples
#'
#' source(system.file("extdata/sample_data.R", package = "nhdplusTools"))
#'
#' nhdplus_path(sample_data)
#'
#' staged_nhdplus <- stage_national_data(output_path = tempdir())
#'
#' sample_flines <- readRDS(staged_nhdplus$flowline)
#'
#' subset_rpu(sample_flines, rpu = "07b")
subset_rpu <- function(fline, rpu, run_make_standalone = TRUE) {
  # Find all outlets of current rpu and sort by size
  # !ToNode %in% FromNode finds non-terminal flowlines that exit the domain.
  outlets <- filter(fline, .data$RPUID %in% rpu)

  if("tocomid" %in% names(outlets)) outlets <- dplyr::rename(outlets, toCOMID = .data$tocomid)

  if(any(c("tocomid", "toCOMID") %in% names(outlets))) {
    outlets <- st_sf(filter(outlets, .data$Hydroseq == .data$TerminalPa |
                              !.data$toCOMID %in% .data$COMID))
  } else {

    outlets <- st_sf(filter(outlets, .data$TerminalFl == 1 |
                              !.data$ToNode %in% .data$FromNode))
  }
  outlets <- arrange(outlets, desc(.data$ArbolateSu))

  # run nhdplusTools::get_UT for all outlets and concatenate.
  network <- lapply(outlets$COMID,
                    function(x, fline) get_UT(fline, x),
                    fline = fline)
  network <- do.call(c, network)

  # Filter so only navigable flowlines are included.
  fline <- fline[fline$COMID %in% network, ]

  # For flowlines labaled as in the RPU, find the top and bottom of each
  # LevelPath. This was required for some unique network situations.
  fline_sub <- filter(drop_geometry(fline), .data$RPUID %in% rpu)

  fline_sub <- group_by(fline_sub, .data$LevelPathI)

  fline_sub <- summarize(fline_sub,
                         lp_top = max(.data$Hydroseq),
                         lp_bot = min(.data$Hydroseq))

  # Using the levelpath top and bottoms found above, filter the complete
  # domain to the hydrosequence of the levelpath top and bottoms instead
  # of trusting the RPUID to be useable.
  fline <- left_join(fline, fline_sub, by = "LevelPathI")

  fline <- group_by(filter(fline, .data$LevelPathI %in% fline_sub$LevelPathI),
                    .data$LevelPathI)

  fline <- ungroup(filter(fline, .data$Hydroseq >= .data$lp_bot &
                            .data$Hydroseq <= .data$lp_top))

  if(run_make_standalone) {
    make_standalone(fline)
  } else {
    fline
  }
}

#' @noRd
get_nhdplus_byid <- function(comids, layer, streamorder = NULL) {

  if(layer == "nhdflowline_network"){
    query_usgs_geoserver(ids = comids, type = "nhd", filter = streamorder_filter(streamorder))
  } else if(layer == "catchmentsp"){
    query_usgs_geoserver(ids = comids, type = "catchment")
  } else {
    stop("Layer must be one of catchmentsp, nhdflowline_network")
  }
}

#' @noRd
get_nhdplus_bybox <- function(box, layer, streamorder = NULL) {

  if(!layer %in% c("nhdarea", "nhdwaterbody", "nhdflowline_network",
                   "nhdflowline_nonnetwork", "catchmentsp")) {
    stop("Layer must be one of nhdarea, nhdwaterbody")
  }

  type <- dplyr::filter(query_usgs_geoserver(),
                        .data$geoserver == layer)$user_call

  query_usgs_geoserver(AOI = box,
                       type = type)
}
