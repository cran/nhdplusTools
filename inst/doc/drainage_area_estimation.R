## ----setup, include = FALSE---------------------------------------------------
library(nhdplusTools)

#pkgdown build uses this
build_local <- (Sys.getenv("BUILD_VIGNETTES") == "TRUE")

#cran package uses this
cran_build <- (Sys.getenv("BUILD_VIGNETTES_CRAN") == "TRUE")

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 12,
  fig.height = 7,
  dpi = 150,
  out.width = "100%",
  eval = build_local && !cran_build,
  fig.path = "drainage_area_figures/"
)

huc12 <- NULL
wbd_gdb <- "../../Data/final_wbd/WBD_National_GDB/WBD_National_GDB.gdb"
if(file.exists(wbd_gdb) & !cran_build) {
  huc12 <- sf::read_sf(wbd_gdb, "WBDHU12")
  # Rename columns to the schema expected by
  # get_drainage_area_estimates() and the vignette plot helpers.
  rename_map <- c(
    huc12 = "huc_12",
    noncontributingareaacres = "ncontrb_a",
    hutype = "hu_12_type"
  )
  for(src in names(rename_map)) {
    if(src %in% names(huc12) && !(rename_map[[src]] %in% names(huc12)))
      names(huc12)[names(huc12) == src] <- rename_map[[src]]
  }
}

# Local NHDPlusV2 national snapshot. Gated behind USE_LOCAL_NHDP so the
# read only happens when explicitly requested. When enabled and the GDB
# is present, NHDWaterbody and CatchmentSP are passed through
# waterbody_data / catchment_data so the point-in-waterbody lookup
# (Malheur), the gap-zone catchment retrieval, and the full-network
# catchment retrieval are resolved locally instead of via the OGC API.
# NHDFlowline_Network is read for the vignette's flowline-plot fetch.
waterbody_data <- NULL
catchment_data <- NULL
flowlines_data <- NULL
if(Sys.getenv("USE_LOCAL_NHDP") == "TRUE") {
  nhdp_gdb <- paste0("../../Data/nhdp/NHDPlusNationalData/",
    "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb")
  if(file.exists(nhdp_gdb)) {
    waterbody_data <- sf::read_sf(nhdp_gdb, "NHDWaterbody")
    catchment_data <- sf::read_sf(nhdp_gdb, "CatchmentSP")
    flowlines_data <- sf::read_sf(nhdp_gdb, "NHDFlowline_Network")
    flowlines_data <- rename_geometry(flowlines_data, "geom")
    names(flowlines_data) <- tolower(names(flowlines_data))
  }
}

has_davidson <- FALSE
oldoption <- options(scipen = 9999)

## ----fetch_data, message = TRUE-----------------------------------------------
# library(sf)
# 
# data_dir <- nhdplusTools_data_dir()
# dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
# 
# # On-network HU12 pour points from the Mainstem Rivers data release
# # (Blodgett 2022, doi:10.5066/P13LNDDQ). Downloaded once and cached so
# # subsequent builds reuse the local file instead of hitting the NLDI.
# outlets_path <- file.path(data_dir, "finalwbd_outlets.gpkg")
# outlets_url <- paste0(
#   "https://prod-is-usgs-sb-prod-publish.s3.amazonaws.com/",
#   "65cbc0b3d34ef4b119cb37e9/finalwbd_outlets.gpkg")
# if(!file.exists(outlets_path)) {
#   message("Downloading HU12 outlets GPKG to ", outlets_path)
#   download.file(outlets_url, outlets_path, mode = "wb")
# }
# 
# # Define the six study sites
# sites <- list(
#   black_earth = list(featureSource = "nwissite",
#     featureID = "USGS-05406500"),
#   french_broad = list(featureSource = "nwissite",
#     featureID = "USGS-03451500"),
#   brazos = list(featureSource = "nwissite",
#     featureID = "USGS-08116650"),
#   james = list(featureSource = "nwissite",
#     featureID = "USGS-06468000"),
#   malheur = st_sfc(st_point(c(-118.8, 43.3)), crs = 4326),
#   purgatoire = list(featureSource = "nwissite",
#     featureID = "USGS-07128500"),
#   davidson = list(featureSource = "nwissite",
#     featureID = "USGS-08110075")
# )
# 
# # Skip NHDPlusHR for large basins to keep fetch times reasonable
# hr_basins <- c("black_earth", "french_broad", "davidson", "purgatoire",
#   "james")
# 
# fetch_or_load <- function(name, start) {
#   rds_path <- file.path(data_dir,
#     paste0("da_est_", name, ".rds"))
#   if(file.exists(rds_path)) {
#     message("Loading cached: ", name)
#     return(readRDS(rds_path))
#   }
#   message("Fetching: ", name)
#   result <- tryCatch(
#     get_drainage_area_estimates(start, catchments = TRUE,
#       huc12_data = huc12,
#       huc12_outlets = outlets_path,
#       waterbody_data = waterbody_data,
#       catchment_data = catchment_data,
#       local_navigation = TRUE,
#       nhdplushr = name %in% hr_basins),
#     error = function(e) {
#       warning("Failed for ", name, ": ", conditionMessage(e),
#         call. = FALSE)
#       NULL
#     })
#   if(!is.null(result)) saveRDS(result, rds_path)
#   result
# }
# 
# da_results <- Map(fetch_or_load, names(sites), sites)
# da_results <- Filter(Negate(is.null), da_results)

## ----fetch_flowlines, message = TRUE------------------------------------------
# fetch_flowlines <- function(name, da_result) {
#   rds_path <- file.path(data_dir, paste0("fl_", name, ".rds"))
#   if(file.exists(rds_path)) {
#     message("Loading cached flowlines: ", name)
#     return(readRDS(rds_path))
#   }
#   message("Fetching flowlines: ", name)
#   comids <- da_result$all_network$comid
#   fl <- if(!is.null(flowlines_data)) {
#     message("  Subsetting from local NHDFlowline_Network...")
#     flowlines_data[flowlines_data$comid %in% comids, ]
#   } else {
#     tryCatch(
#       get_nhdplus(comid = comids, realization = "flowline"),
#       error = function(e) {
#         warning("Flowline fetch failed for ", name, ": ",
#           conditionMessage(e), call. = FALSE)
#         NULL
#       })
#   }
#   if(!is.null(fl) && nrow(fl) > 0) saveRDS(fl, rds_path)
#   fl
# }
# 
# flowlines <- Map(fetch_flowlines, names(da_results), da_results)
# flowlines <- Filter(Negate(is.null), flowlines)

## ----fetch_nwis_da, message = TRUE--------------------------------------------
# # Fetch NWIS drainage area for sites with an nwissite featureSource
# nwis_ids <- vapply(sites, function(s) {
#   if(is.list(s) && identical(s$featureSource, "nwissite"))
#     s$featureID else NA_character_
# }, character(1))
# nwis_ids <- nwis_ids[!is.na(nwis_ids)]
# 
# nwis_da <- if(length(nwis_ids) > 0) {
#   tryCatch({
#     ml <- dataRetrieval::read_waterdata_monitoring_location(
#       unname(nwis_ids))
#     # Match returned rows back to site names by monitoring_location_id
#     idx <- match(nwis_ids, ml$monitoring_location_id)
#     # drainage_area and contributing_drainage_area are in sq miles
#     data.frame(
#       name = names(nwis_ids),
#       nwis_da_sqmi = ml$drainage_area[idx],
#       nwis_contrib_da_sqmi = ml$contributing_drainage_area[idx],
#       nwis_da_sqkm = ml$drainage_area[idx] * 2.58999,
#       nwis_contrib_da_sqkm = ml$contributing_drainage_area[idx] *
#         2.58999,
#       stringsAsFactors = FALSE
#     )
#   }, error = function(e) {
#     warning("NWIS site fetch failed: ", conditionMessage(e),
#       call. = FALSE)
#     NULL
#   })
# } else NULL

## ----structure----------------------------------------------------------------
# names(da_results$black_earth)

## ----summary_table------------------------------------------------------------
# basin_labels <- c(
#   black_earth = "Black Earth Creek",
#   french_broad = "French Broad",
#   brazos = "Brazos at Rosharon",
#   james = "James River",
#   malheur = "Malheur Lake",
#   davidson = "Davidson Creek",
#   purgatoire = "Purgatoire River"
# )
# 
# summary_df <- data.frame(
#   basin = basin_labels[names(da_results)],
#   network_da = vapply(da_results, \(x) x$network_da_sqkm, numeric(1)),
#   da_huc12 = vapply(da_results, \(x) x$da_huc12_sqkm, numeric(1)),
#   da_huc10 = vapply(da_results,
#     \(x) ifelse(is.na(x$da_huc10_sqkm), NA_real_, x$da_huc10_sqkm),
#     numeric(1)),
#   da_huc08 = vapply(da_results,
#     \(x) ifelse(is.na(x$da_huc08_sqkm), NA_real_, x$da_huc08_sqkm),
#     numeric(1)),
#   nhdplushr = vapply(da_results,
#     \(x) ifelse(is.na(x$nhdplushr_network_dasqkm), NA_real_,
#       x$nhdplushr_network_dasqkm),
#     numeric(1))
# )
# 
# # Add NWIS drainage areas where available
# if(!is.null(nwis_da)) {
#   summary_df$nwis_da <- ifelse(
#     names(da_results) %in% nwis_da$name,
#     nwis_da$nwis_da_sqkm[match(names(da_results), nwis_da$name)],
#     NA_real_)
#   summary_df$nwis_contrib_da <- ifelse(
#     names(da_results) %in% nwis_da$name,
#     nwis_da$nwis_contrib_da_sqkm[match(names(da_results), nwis_da$name)],
#     NA_real_)
# } else {
#   summary_df$nwis_da <- NA_real_
#   summary_df$nwis_contrib_da <- NA_real_
# }
# 
# knitr::kable(summary_df, digits = 1,
#   col.names = c("Basin", "Network DA", "HU12 DA",
#     "HU10 DA", "HU8 DA", "NHDPlusHR DA",
#     "NWIS DA", "NWIS Contributing DA"),
#   caption = "Drainage area estimates (sq km) by basin and method")

## ----helpers, include = FALSE-------------------------------------------------
# library(ggplot2)
# 
# # USGS USTopo provider for maptiles. The "?f=jpg" suffix is a no-op
# # query parameter that gives maptiles a recognized extension hint.
# usgs_topo_provider <- list(
#   src = "USGS_USTopo",
#   q = paste0("https://basemap.nationalmap.gov/arcgis/rest/services/",
#     "USGSTopo/MapServer/tile/{z}/{y}/{x}?f=jpg"),
#   sub = NA,
#   cit = "USGS The National Map: USGSTopo"
# )
# 
# # Cache directory for tiles
# tile_cache <- file.path(nhdplusTools_data_dir(), "da_v_tiles")
# dir.create(tile_cache, recursive = TRUE, showWarnings = FALSE)
# 
# # Compute a target tile zoom level so the rendered tiles are sharp
# # at the figure size we use. Web Mercator world width is 40075016 m
# # and tile size is 256 px. The figure pixel width drives the
# # resolution we want from the tile mosaic.
# target_zoom <- function(sf_obj, fig_px = 1700) {
#   bb <- st_bbox(sf_obj)
#   span_m <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
#   if(!is.finite(span_m) || span_m <= 0) return(10L)
#   z <- log2(40075016 * fig_px / (256 * span_m))
#   as.integer(min(max(round(z), 4L), 15L))
# }
# 
# # Fetch USGS topo tiles for an sf bbox. Returns a list with a color
# # matrix (suitable for annotation_raster) and the tile bounding box.
# # Uses annotation_raster (not geom_raster) so the fill aesthetic
# # stays free for choropleth use. When fig_ratio is provided, the
# # bbox is expanded to match the figure aspect ratio via
# # fit_bbox_to_fig so that tiles cover the full canvas.
# fetch_topo_raster <- function(sf_obj, zoom = NULL,
#   fig_ratio = NULL) {
#   fit <- if(!is.null(fig_ratio))
#     fit_bbox_to_fig(sf_obj, fig_ratio = fig_ratio) else NULL
#   if(!is.null(fit)) {
#     bb <- st_bbox(c(
#       xmin = fit$xlim[1], xmax = fit$xlim[2],
#       ymin = fit$ylim[1], ymax = fit$ylim[2]
#     ), crs = st_crs(3857))
#     bb_sf <- st_as_sfc(bb)
#   } else {
#     bb <- st_bbox(sf_obj)
#     pad_x <- (bb["xmax"] - bb["xmin"]) * 0.05
#     pad_y <- (bb["ymax"] - bb["ymin"]) * 0.05
#     bb["xmin"] <- bb["xmin"] - pad_x
#     bb["xmax"] <- bb["xmax"] + pad_x
#     bb["ymin"] <- bb["ymin"] - pad_y
#     bb["ymax"] <- bb["ymax"] + pad_y
#     bb_sf <- st_as_sfc(bb, crs = st_crs(sf_obj))
#   }
#   if(is.null(zoom))
#     zoom <- target_zoom(st_transform(bb_sf, 3857))
#   args <- list(
#     x = bb_sf, provider = usgs_topo_provider, crop = TRUE,
#     cachedir = tile_cache, zoom = zoom
#   )
#   tiles <- tryCatch(
#     do.call(maptiles::get_tiles, args),
#     error = function(e) {
#       warning("Tile fetch failed: ", conditionMessage(e),
#         call. = FALSE)
#       NULL
#     })
#   if(is.null(tiles)) return(NULL)
#   rgb_arr <- terra::as.array(tiles)
#   if(length(dim(rgb_arr)) != 3 || dim(rgb_arr)[3] < 3) return(NULL)
#   col_mat <- matrix(
#     rgb(rgb_arr[, , 1], rgb_arr[, , 2], rgb_arr[, , 3],
#       maxColorValue = 255),
#     nrow = dim(rgb_arr)[1], ncol = dim(rgb_arr)[2]
#   )
#   ext <- terra::ext(tiles)
#   list(
#     raster = col_mat,
#     xmin = ext[1], xmax = ext[2],
#     ymin = ext[3], ymax = ext[4]
#   )
# }
# 
# # Expand a bbox to match the figure aspect ratio so that coord_sf
# # fills the canvas without white sidebars. The basin content is
# # centered and the narrower axis is padded with additional basemap.
# fit_bbox_to_fig <- function(sf_obj, fig_ratio = 12 / 7, pad = 0.08) {
#   bb <- st_bbox(st_transform(sf_obj, 3857))
#   dx <- unname(bb["xmax"] - bb["xmin"])
#   dy <- unname(bb["ymax"] - bb["ymin"])
#   if(!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0)
#     return(NULL)
#   # Add uniform padding first
#   cx <- unname((bb["xmin"] + bb["xmax"]) / 2)
#   cy <- unname((bb["ymin"] + bb["ymax"]) / 2)
#   dx <- dx * (1 + pad * 2)
#   dy <- dy * (1 + pad * 2)
#   # Expand the narrower axis to match the target aspect ratio
#   current_ratio <- dx / dy
#   if(current_ratio < fig_ratio) {
#     dx <- dy * fig_ratio
#   } else {
#     dy <- dx / fig_ratio
#   }
#   list(
#     xlim = c(cx - dx / 2, cx + dx / 2),
#     ylim = c(cy - dy / 2, cy + dy / 2)
#   )
# }
# 
# # Add basemap as annotation_raster (does not consume fill scale)
# add_topo <- function(p, sf_obj, fig_ratio = 12 / 7) {
#   topo <- fetch_topo_raster(sf_obj, fig_ratio = fig_ratio)
#   if(is.null(topo)) return(p)
#   p + annotation_raster(
#     topo$raster,
#     xmin = topo$xmin, xmax = topo$xmax,
#     ymin = topo$ymin, ymax = topo$ymax
#   )
# }
# 
# # Resolve gage point in target CRS
# get_gage_pt <- function(da, target_crs) {
#   pt <- st_geometry(da$start_feature)
#   if(!inherits(pt, "sfc_POINT"))
#     pt <- st_centroid(pt)
#   st_transform(st_sf(geometry = pt), target_crs)
# }
# 
# # Pick the spatially-largest of the available HUC12 sets returned by
# # get_drainage_area_estimates: hu12_by_huc08 (if present), otherwise
# # hu12_by_huc10, otherwise hu12_by_huc12. Largest is determined by
# # the sum of HUC12 dasqkm in each set.
# largest_hu12_set <- function(da) {
#   candidates <- list(
#     huc12 = da$hu12_by_huc12,
#     huc10 = da$hu12_by_huc10,
#     huc08 = da$hu12_by_huc08
#   )
#   candidates <- Filter(
#     function(x) !is.null(x) && nrow(x) > 0, candidates)
#   areas <- vapply(candidates,
#     function(x) sum(x$dasqkm, na.rm = TRUE), numeric(1))
#   candidates[[which.max(areas)]]
# }
# 
# # Compute a single merged basin polygon for a given HUC level set
# # (hu12_by_huc12, hu12_by_huc10, or hu12_by_huc08), unioned with
# # extra_catchments and the local catchments around the gage.
# # Uses hydroloom::dissolve_polygons() with gap_tolerance to absorb
# # sliver gaps between polygons from different services (WBD HUC12s,
# # NHDPlusV2 catchments, split-catchment service).
# compute_basin_polygon <- function(da, hu_layer, target_crs,
#   tol_m = 250) {
#   common_crs <- st_crs(hu_layer)
#   pieces <- list()
#   pieces$hu <- st_geometry(hu_layer)
# 
#   if(!is.null(da$extra_catchments) && nrow(da$extra_catchments) > 0)
#     pieces$extra <- st_geometry(
#       st_transform(da$extra_catchments, common_crs))
# 
#   if(!is.null(da$split_catchment) && nrow(da$split_catchment) > 0) {
#     sc <- da$split_catchment
#     sc_full <- sc[sc$id == "catchment", ]
#     if(nrow(sc_full) > 0)
#       pieces$catch <- st_geometry(st_transform(sc_full, common_crs))
#   }
# 
#   all_polys <- st_sf(
#     geometry = do.call(c, pieces)
#   ) |> st_set_crs(common_crs)
# 
#   hydroloom::dissolve_polygons(all_polys,
#     gap_tolerance = tol_m,
#     max_hole_area = Inf,
#     single_polygon = TRUE) |>
#     st_transform(target_crs)
# }
# 
# # Common theme: no axis ticks, labels, or grid; legend sized for
# # readability
# map_theme <- function() {
#   theme_void() +
#     theme(
#       legend.position = "bottom",
#       legend.box = "horizontal",
#       legend.text = element_text(size = 11),
#       legend.title = element_text(size = 12),
#       legend.key.size = unit(0.6, "cm"),
#       plot.title = element_blank(),
#       plot.margin = margin(2, 2, 2, 2)
#     )
# }
# 
# # Add the gray basin underlay to a plot. The polygon represents the
# # full estimated drainage area envelope (largest available HUC level
# # plus extra catchments and local catchments).
# add_basin_underlay <- function(p, basin_poly) {
#   p + geom_sf(data = basin_poly,
#     fill = "gray30", color = NA, alpha = 0.15,
#     inherit.aes = FALSE)
# }
# 
# # Mode 1: Drainage Area Boundaries
# plot_boundaries <- function(da, basin_label) {
#   # Reproject working layers to Web Mercator (3857) so the basemap
#   # tiles align cleanly.
#   target_crs <- st_crs(3857)
#   gage_pt <- get_gage_pt(da, target_crs)
#   basin_poly <- compute_basin_polygon(da, largest_hu12_set(da), target_crs)
# 
#   # Dissolve boundaries at each HU level. Each HU level is combined
#   # with extra_catchments and split_catchment via compute_basin_polygon
#   # so every boundary reaches the gage instead of stopping at the
#   # lowest complete HU outlet.
#   boundaries <- list()
#   boundaries[["HU12"]] <- st_geometry(
#     compute_basin_polygon(da, da$hu12_by_huc12, target_crs))
#   if(!is.null(da$hu12_by_huc10))
#     boundaries[["HU10"]] <- st_geometry(
#       compute_basin_polygon(da, da$hu12_by_huc10, target_crs))
#   if(!is.null(da$hu12_by_huc08))
#     boundaries[["HU8"]] <- st_geometry(
#       compute_basin_polygon(da, da$hu12_by_huc08, target_crs))
#   if(!is.null(da$all_catchments))
#     boundaries[["Network catchments (NHDPlusV2)"]] <- st_union(
#       st_transform(da$all_catchments, target_crs))
#   if(!is.null(da$nhdplushr_boundary))
#     boundaries[["NHDPlusHR"]] <- st_transform(
#       da$nhdplushr_boundary, target_crs)
# 
#   # Linetype, color, and linewidth all keyed off the source factor so
#   # one merged "DA boundary" legend distinguishes the layers.
#   # Network catchments (NHDPlusV2) gets a saturated magenta dotted
#   # line so it reads distinctly against the basemap and against HU10's
#   # longdash. NHDPlusHR uses orange to separate from the HU layers.
#   all_ltypes <- c(
#     "HU12" = "solid",
#     "HU10" = "longdash",
#     "HU8" = "dotdash",
#     "Network catchments (NHDPlusV2)" = "22",
#     "NHDPlusHR" = "solid"
#   )
#   all_colors <- c(
#     "HU12" = "gray10",
#     "HU10" = "gray10",
#     "HU8" = "gray10",
#     "Network catchments (NHDPlusV2)" = "magenta",
#     "NHDPlusHR" = "#D55E00"
#   )
#   all_lwds <- c(
#     "HU12" = 0.8,
#     "HU10" = 0.8,
#     "HU8" = 0.8,
#     "Network catchments (NHDPlusV2)" = 0.5,
#     "NHDPlusHR" = 0.6
#   )
#   ltypes  <- all_ltypes[names(boundaries)]
#   lcolors <- all_colors[names(boundaries)]
#   lwds    <- all_lwds[names(boundaries)]
# 
#   bnd_sf <- do.call(rbind, lapply(seq_along(boundaries), function(i) {
#     st_sf(source = factor(names(boundaries)[i],
#         levels = names(boundaries)),
#       geometry = boundaries[[i]])
#   }))
# 
#   hu12_pp <- if(!is.null(da$hu12_outlet))
#     st_transform(da$hu12_outlet, target_crs) else NULL
# 
#   # Combine gage + HU12 outlets into one points layer with a single
#   # shape scale. Keeping shape only (no color aes) frees the color
#   # scale for the boundary symbology.
#   pts_list <- list(
#     st_sf(label = "Stream gage / outlet",
#       geometry = st_geometry(gage_pt))
#   )
#   if(!is.null(hu12_pp))
#     pts_list <- c(pts_list, list(
#       st_sf(label = "HU12 outlet",
#         geometry = st_geometry(hu12_pp))))
#   pts_sf <- do.call(rbind, pts_list)
#   pts_sf$label <- factor(pts_sf$label,
#     levels = c("Stream gage / outlet", "HU12 outlet"))
# 
#   fig_lims <- fit_bbox_to_fig(basin_poly)
# 
#   p <- ggplot() |>
#     add_topo(basin_poly) |>
#     add_basin_underlay(basin_poly)
# 
#   # HU12 outlets drawn first so basin boundaries overlay them.
#   if(!is.null(hu12_pp))
#     p <- p + geom_sf(
#       data = pts_sf[pts_sf$label == "HU12 outlet", ],
#       shape = 4, size = 1, color = "darkred", stroke = 0.6,
#       alpha = 0.75, inherit.aes = FALSE)
# 
#   p <- p +
#     geom_sf(data = bnd_sf,
#       aes(linetype = source, color = source, linewidth = source),
#       fill = NA, inherit.aes = FALSE) +
#     scale_linetype_manual(values = ltypes,  name = "DA boundary") +
#     scale_color_manual(   values = lcolors, name = "DA boundary") +
#     scale_linewidth_manual(values = lwds,   name = "DA boundary") +
#     geom_sf(
#       data = pts_sf[pts_sf$label == "Stream gage / outlet", ],
#       shape = 17, size = 5, color = "black", fill = "white",
#       stroke = 1.4, inherit.aes = FALSE) +
#     coord_sf(crs = target_crs, expand = FALSE,
#       xlim = fig_lims$xlim, ylim = fig_lims$ylim) +
#     labs(title = paste(basin_label, "-- Drainage Area Boundaries")) +
#     map_theme()
# 
#   p
# }
# 
# # Mode 2: Stream Network
# plot_network <- function(da, fl, basin_label) {
#   target_crs <- st_crs(3857)
#   gage_pt <- get_gage_pt(da, target_crs)
#   fl <- st_transform(fl, target_crs)
#   basin_poly <- compute_basin_polygon(da, largest_hu12_set(da), target_crs)
# 
#   # fcode 46003 = intermittent, 46007 = ephemeral, 46006 = perennial
#   fl$flow_class <- ifelse(
#     fl$fcode %in% c(46003L, 46007L), "Intermittent/Ephemeral",
#     "Perennial/Other"
#   )
#   fl$lwd <- pmin(pmax(fl$streamorde, 1) / 3, 1.5)
# 
#   pts <- st_sf(label = "Stream gage / outlet",
#     geometry = st_geometry(gage_pt))
# 
#   fig_lims <- fit_bbox_to_fig(basin_poly)
# 
#   ggplot() |>
#     add_topo(basin_poly) |>
#     add_basin_underlay(basin_poly) +
#     geom_sf(data = fl,
#       aes(color = flow_class, linewidth = lwd),
#       show.legend = "line", inherit.aes = FALSE) +
#     scale_color_manual(
#       values = c("Perennial/Other" = "steelblue",
#         "Intermittent/Ephemeral" = "lightskyblue"),
#       name = "Flow class") +
#     scale_linewidth_identity() +
#     ggnewscale_color(pts, target_crs) +
#     coord_sf(crs = target_crs, expand = FALSE,
#       xlim = fig_lims$xlim, ylim = fig_lims$ylim) +
#     labs(title = paste(basin_label, "-- Stream Network")) +
#     map_theme()
# }
# 
# # Brazos-specific variant of plot_network that adds a Caprock Escarpment
# # trace and label per LH-217625600 / LH-2005481448. The escarpment is
# # the eastern margin of the Llano Estacado and the physiographic
# # boundary between the High Plains and the Rolling Plains (Wermund,
# # 1996). The polyline below approximates its trace through the Brazos
# # basin's northwestern arm.
# plot_network_brazos <- function(da, fl, basin_label) {
#   target_crs <- st_crs(3857)
#   caprock_coords <- matrix(c(
#     -101.05, 34.50,
#     -101.00, 33.70,
#     -100.95, 32.85
#   ), ncol = 2, byrow = TRUE)
#   caprock_line <- sf::st_transform(
#     sf::st_sfc(sf::st_linestring(caprock_coords), crs = 4326),
#     target_crs)
#   caprock_anchor <- sf::st_transform(
#     sf::st_sfc(sf::st_point(c(-101.05, 34.50)), crs = 4326),
#     target_crs)
#   caprock_line_sf <- sf::st_sf(geometry = caprock_line)
#   caprock_label_sf <- sf::st_sf(
#     label = "Caprock Escarpment",
#     geometry = caprock_anchor)
# 
#   plot_network(da, fl, basin_label) +
#     geom_sf(data = caprock_line_sf,
#       color = "black", linewidth = 0.8,
#       inherit.aes = FALSE) +
#     geom_sf_text(data = caprock_label_sf, aes(label = label),
#       size = 3.8, fontface = "italic", color = "gray15",
#       hjust = 0, nudge_x = 8000, nudge_y = 22000,
#       inherit.aes = FALSE)
# }
# 
# # Helper to add the gage point as a labeled layer using a shape scale
# # that doesn't conflict with other color scales already in use.
# ggnewscale_color <- function(pts, target_crs) {
#   list(
#     geom_sf(data = pts, aes(shape = label),
#       size = 5, color = "black", fill = "white", stroke = 1.4,
#       inherit.aes = FALSE),
#     scale_shape_manual(
#       values = c("Stream gage / outlet" = 17),
#       name = NULL)
#   )
# }
# 
# # Mode 3: HUC12 Type (hu_12_type)
# plot_type <- function(da, fl, basin_label,
#   fig_ratio = 12 / 8) {
#   target_crs <- st_crs(3857)
#   largest <- largest_hu12_set(da)
#   hu12 <- st_transform(largest, target_crs)
# 
#   type_labels <- c(
#     "S" = "Standard",
#     "C" = "Closed basin",
#     "M" = "Multiple outlets",
#     "F" = "Frontal",
#     "W" = "Water",
#     "D" = "Drainage",
#     "I" = "Island"
#   )
#   hu12$type_label <- factor(
#     type_labels[hu12$hu_12_type],
#     levels = c("Standard", "Closed basin", "Multiple outlets",
#       "Frontal", "Water", "Drainage", "Island")
#   )
# 
#   gage_pt <- get_gage_pt(da, target_crs)
#   fl <- st_transform(fl, target_crs)
#   basin_poly <- compute_basin_polygon(da, largest, target_crs)
# 
#   pts <- st_sf(label = "Stream gage / outlet",
#     geometry = st_geometry(gage_pt))
# 
#   fig_lims <- fit_bbox_to_fig(basin_poly, fig_ratio = fig_ratio)
# 
#   # 30% alpha baked into the fill so the legend swatches match the map.
#   # Standard is light grey so the legend swatch reads against white.
#   col30 <- function(col) adjustcolor(col, alpha.f = 0.3)
#   type_colors <- c(
#     "Standard" = adjustcolor("gray85", alpha.f = 0.4),
#     "Closed basin" = col30("#8B4513"),
#     "Multiple outlets" = col30("#555555"),
#     "Frontal" = col30("#228B22"),
#     "Water" = col30("#1E90FF"),
#     "Drainage" = col30("#9467BD"),
#     "Island" = col30("#FFC107")
#   )
# 
#   p <- ggplot() |>
#     add_topo(basin_poly, fig_ratio = fig_ratio) |>
#     add_basin_underlay(basin_poly) +
#     geom_sf(data = fl, color = "steelblue", linewidth = 0.25,
#       alpha = 0.5, inherit.aes = FALSE) +
#     geom_sf(data = hu12,
#       aes(fill = type_label),
#       color = "gray50", linewidth = 0.5,
#       inherit.aes = FALSE) +
#     scale_fill_manual(
#       values = type_colors,
#       drop = TRUE,
#       name = "HU12 type")
# 
#   p +
#     ggnewscale_color(pts, target_crs) +
#     coord_sf(crs = target_crs, expand = FALSE,
#       xlim = fig_lims$xlim, ylim = fig_lims$ylim) +
#     labs(title = paste(basin_label, "-- HU12 Types")) +
#     map_theme() +
#     theme(legend.box = "vertical") +
#     guides(
#       fill = guide_legend(
#         nrow = 1, order = 1,
#         override.aes = list(color = "gray50", linewidth = 0.5)))
# }
# 
# # Per-basin summary table
# basin_summary_table <- function(da, basin_name = NULL) {
#   vals <- data.frame(
#     source = c("Network", "HU12",
#       "HU10", "HU8", "NHDPlusHR"),
#     total_sqkm = c(
#       da$network_da_sqkm,
#       da$da_huc12_sqkm,
#       ifelse(is.na(da$da_huc10_sqkm), NA_real_, da$da_huc10_sqkm),
#       ifelse(is.na(da$da_huc08_sqkm), NA_real_, da$da_huc08_sqkm),
#       ifelse(is.na(da$nhdplushr_network_dasqkm), NA_real_,
#         da$nhdplushr_network_dasqkm)
#     )
#   )
#   # Add NWIS rows if this basin has NWIS data
#   if(!is.null(nwis_da) && !is.null(basin_name) &&
#     basin_name %in% nwis_da$name) {
#     row <- nwis_da[nwis_da$name == basin_name, ]
#     vals <- rbind(vals, data.frame(
#       source = c("NWIS", "NWIS contributing"),
#       total_sqkm = c(row$nwis_da_sqkm, row$nwis_contrib_da_sqkm)
#     ))
#   }
#   knitr::kable(vals, digits = 1,
#     col.names = c("Source", "Area (sq km)"),
#     caption = "Drainage area estimates by source")
# }

## ----french_broad_boundaries--------------------------------------------------
# plot_boundaries(da_results$french_broad, "French Broad")

## ----fb_table-----------------------------------------------------------------
# basin_summary_table(da_results$french_broad, "french_broad")

## ----french_broad_network-----------------------------------------------------
# plot_network(da_results$french_broad, flowlines$french_broad,
#   "French Broad")

## ----french_broad_type, fig.height = 8----------------------------------------
# plot_type(da_results$french_broad, flowlines$french_broad,
#   "French Broad")

## ----brazos_boundaries--------------------------------------------------------
# plot_boundaries(da_results$brazos, "Brazos at Rosharon")

## ----bz_table-----------------------------------------------------------------
# basin_summary_table(da_results$brazos, "brazos")

## ----brazos_network-----------------------------------------------------------
# plot_network_brazos(da_results$brazos, flowlines$brazos,
#   "Brazos at Rosharon")

## ----brazos_type, fig.height = 8----------------------------------------------
# plot_type(da_results$brazos, flowlines$brazos,
#   "Brazos at Rosharon")

## ----james_boundaries---------------------------------------------------------
# plot_boundaries(da_results$james, "James River")

## ----jm_table-----------------------------------------------------------------
# basin_summary_table(da_results$james, "james")

## ----james_network------------------------------------------------------------
# plot_network(da_results$james, flowlines$james, "James River")

## ----james_type, fig.height = 8-----------------------------------------------
# plot_type(da_results$james, flowlines$james, "James River")

## ----black_earth_boundaries---------------------------------------------------
# plot_boundaries(da_results$black_earth, "Black Earth Creek")

## ----be_table-----------------------------------------------------------------
# basin_summary_table(da_results$black_earth, "black_earth")

## ----black_earth_network------------------------------------------------------
# plot_network(da_results$black_earth, flowlines$black_earth,
#   "Black Earth Creek")

## ----black_earth_type, fig.height = 8-----------------------------------------
# plot_type(da_results$black_earth, flowlines$black_earth,
#   "Black Earth Creek")

## ----malheur_boundaries-------------------------------------------------------
# plot_boundaries(da_results$malheur, "Malheur Lake")

## ----ml_table-----------------------------------------------------------------
# basin_summary_table(da_results$malheur, "malheur")

## ----malheur_network----------------------------------------------------------
# plot_network(da_results$malheur, flowlines$malheur, "Malheur Lake")

## ----malheur_type, fig.height = 8---------------------------------------------
# plot_type(da_results$malheur, flowlines$malheur, "Malheur Lake")

## ----purgatoire_boundaries----------------------------------------------------
# plot_boundaries(da_results$purgatoire, "Purgatoire River")

## ----purgatoire_table---------------------------------------------------------
# basin_summary_table(da_results$purgatoire, "purgatoire")

## ----purgatoire_network-------------------------------------------------------
# plot_network(da_results$purgatoire, flowlines$purgatoire,
#   "Purgatoire River")

## ----purgatoire_type, fig.height = 8------------------------------------------
# plot_type(da_results$purgatoire, flowlines$purgatoire,
#   "Purgatoire River")

## ----outlet_split_prep, include = FALSE---------------------------------------
# da_dav <- da_results$davidson
# fl_dav <- flowlines$davidson
# has_davidson <- !is.null(da_dav) && !is.null(fl_dav) &&
#   !is.null(da_dav$outlet_split_catchment)
# 
# if(has_davidson) {
#   target_crs <- st_crs(3857)
# 
#   gage_pt <- get_gage_pt(da_dav, target_crs)
#   gage_comid <- as.integer(da_dav$start_feature$comid)
# 
#   # all catchments in the gap between the gage and HUC12 outlets
#   extra_cat <- st_transform(da_dav$extra_catchments, target_crs)
#   all_cat <- st_transform(da_dav$all_catchments, target_crs)
#   fl_sf <- st_transform(fl_dav, target_crs)
#   hu12_pts <- st_transform(da_dav$hu12_outlet, target_crs)
# 
#   # HUC12 split catchment (split at the HUC12 pour point)
#   hu12_sc <- st_transform(da_dav$split_catchment, target_crs)
#   hu12_catch_full <- hu12_sc[hu12_sc$id == "catchment", ]
#   hu12_catch_split <- hu12_sc[hu12_sc$id == "splitCatchment", ]
#   hu12_comid <- as.integer(na.omit(hu12_sc$catchmentID))
# 
#   # gage outlet split catchment (split at the gage point)
#   osc <- st_transform(da_dav$outlet_split_catchment, target_crs)
#   osc_full <- osc[osc$id == "catchment", ]
#   osc_split <- osc[osc$id == "splitCatchment", ]
# 
#   # flowlines in the gap area
#   gap_comids <- c(gage_comid, extra_cat$featureid, hu12_comid)
#   gap_fl <- fl_sf[fl_sf$comid %in% gap_comids, ]
# 
#   # focus bbox: union of gap catchments + split catchments + gage
#   focus_geom <- st_union(c(
#     st_geometry(extra_cat),
#     st_geometry(hu12_catch_full),
#     st_geometry(osc_full),
#     st_geometry(gage_pt)
#   ))
#   focus_bb <- st_bbox(focus_geom)
#   pad_x <- (focus_bb["xmax"] - focus_bb["xmin"]) * 0.05
#   pad_y <- (focus_bb["ymax"] - focus_bb["ymin"]) * 0.05
#   focus_xlim <- c(focus_bb["xmin"] - pad_x, focus_bb["xmax"] + pad_x)
#   focus_ylim <- c(focus_bb["ymin"] - pad_y, focus_bb["ymax"] + pad_y)
# }

## ----outlet_split_overview, fig.height = 10, eval = has_davidson--------------
# p_overview <- ggplot() |>
#   add_topo(focus_geom) +
#   # all catchments as light underlay
#   geom_sf(data = all_cat, fill = "gray30", color = NA, alpha = 0.12,
#     inherit.aes = FALSE) +
#   # gap catchments
#   geom_sf(data = extra_cat, fill = "lightblue", color = "steelblue",
#     linewidth = 0.3, alpha = 0.5, inherit.aes = FALSE) +
#   # HUC12 split catchment — full outline
#   geom_sf(data = hu12_catch_full, fill = NA, color = "gray10",
#     linewidth = 0.6, linetype = "dashed", inherit.aes = FALSE) +
#   # gage outlet catchment — full outline
#   geom_sf(data = osc_full, fill = NA, color = "gray10",
#     linewidth = 0.6, inherit.aes = FALSE) +
#   # flowlines in the gap
#   geom_sf(data = gap_fl, color = "steelblue", linewidth = 0.5,
#     inherit.aes = FALSE) +
#   # HUC12 outlet
#   geom_sf(data = hu12_pts[hu12_pts$comid == hu12_comid, ],
#     shape = 4, color = "darkred", size = 1, stroke = 0.6, alpha = 0.5,
#     inherit.aes = FALSE) +
#   # gage point
#   geom_sf(data = gage_pt, shape = 17, color = "black",
#     fill = "white", size = 5, stroke = 1.4,
#     inherit.aes = FALSE) +
#   coord_sf(crs = target_crs, xlim = focus_xlim, ylim = focus_ylim,
#     expand = FALSE) +
#   labs(title = paste0("Davidson Creek (USGS streamgage 08110075)",
#     " -- Gage and HU12 Outlet Positions")) +
#   map_theme()
# print(p_overview)

## ----outlet_split_huc_zoom, fig.height = 10, eval = has_davidson--------------
# hu12_outlet_pt <- hu12_pts[hu12_pts$comid == hu12_comid, ]
# hu12_coords <- st_coordinates(hu12_outlet_pt)
# half_side <- 250 # meters in EPSG:3857
# zoom_xlim <- c(hu12_coords[1, "X"] - half_side,
#   hu12_coords[1, "X"] + half_side)
# zoom_ylim <- c(hu12_coords[1, "Y"] - half_side,
#   hu12_coords[1, "Y"] + half_side)
# 
# p_huc_zoom <- ggplot() |>
#   add_topo(hu12_outlet_pt) +
#   # gap catchments
#   geom_sf(data = extra_cat, fill = "lightblue", color = "steelblue",
#     linewidth = 0.3, alpha = 0.5, inherit.aes = FALSE) +
#   # HUC12 split catchment — full outline
#   geom_sf(data = hu12_catch_full, fill = NA, color = "gray10",
#     linewidth = 0.6, linetype = "dashed", inherit.aes = FALSE) +
#   # HUC12 split portion
#   geom_sf(data = hu12_catch_split, fill = "orange", color = "gray10",
#     linewidth = 0.4, alpha = 0.5, inherit.aes = FALSE) +
#   # flowlines in the gap
#   geom_sf(data = gap_fl, color = "steelblue", linewidth = 0.5,
#     inherit.aes = FALSE) +
#   # HUC12 outlet
#   geom_sf(data = hu12_outlet_pt,
#     shape = 4, color = "darkred", size = 1, stroke = 0.6, alpha = 0.5,
#     inherit.aes = FALSE) +
#   coord_sf(crs = target_crs, xlim = zoom_xlim, ylim = zoom_ylim,
#     expand = FALSE) +
#   labs(title = paste0("Davidson Creek -- HU12 Outlet Detail",
#     " (500 m view)")) +
#   map_theme()
# print(p_huc_zoom)

## ----outlet_split_detail, fig.height = 10, eval = has_davidson----------------
# # Build an sf with labeled polygons for a single legend
# split_layers <- rbind(
#   st_sf(
#     role = "Upstream of gage (included)",
#     geometry = st_geometry(osc_split)),
#   st_sf(
#     role = "Downstream of gage (excluded)",
#     geometry = st_difference(
#       st_geometry(osc_full), st_geometry(osc_split))),
#   st_sf(
#     role = "Upstream of HU outlet (excluded, overlaps HU)",
#     geometry = st_geometry(hu12_catch_split)),
#   st_sf(
#     role = "Downstream of HU outlet (included as local area)",
#     geometry = st_difference(
#       st_geometry(hu12_catch_full), st_geometry(hu12_catch_split)))
# )
# 
# split_colors <- c(
#   "Upstream of gage (included)" = "#4DAF4A",
#   "Downstream of gage (excluded)" = "#E41A1C",
#   "Upstream of HU outlet (excluded, overlaps HU)" = "#FF7F00",
#   "Downstream of HU outlet (included as local area)" = "#377EB8"
# )
# 
# split_layers$role <- factor(split_layers$role,
#   levels = names(split_colors))
# 
# p_split <- ggplot() |>
#   add_topo(focus_geom) +
#   # gap catchments as underlay
#   geom_sf(data = extra_cat, fill = "gray80", color = "gray60",
#     linewidth = 0.2, alpha = 0.3, inherit.aes = FALSE) +
#   # split polygons with role-based fill
#   geom_sf(data = split_layers,
#     aes(fill = role), color = "gray20", linewidth = 0.4,
#     alpha = 0.6, inherit.aes = FALSE) +
#   scale_fill_manual(values = split_colors, name = NULL) +
#   # flowlines
#   geom_sf(data = gap_fl, color = "steelblue", linewidth = 0.5,
#     inherit.aes = FALSE) +
#   # HUC12 outlet
#   geom_sf(data = hu12_pts[hu12_pts$comid == hu12_comid, ],
#     shape = 4, color = "darkred", size = 1, stroke = 0.6, alpha = 0.5,
#     inherit.aes = FALSE) +
#   # gage
#   geom_sf(data = gage_pt, shape = 17, color = "black",
#     fill = "white", size = 5, stroke = 1.4,
#     inherit.aes = FALSE) +
#   coord_sf(crs = target_crs, xlim = focus_xlim, ylim = focus_ylim,
#     expand = FALSE) +
#   labs(title = paste0("Davidson Creek -- Split Catchment Roles"),
#     subtitle = paste0(
#       "Flowline measure at gage: ",
#       round(da_dav$outlet_flowline_measure, 1),
#       "; downstream removed: ",
#       round(osc_full$dasqkm - osc_split$dasqkm, 2), " km\u00B2")) +
#   map_theme() +
#   guides(fill = guide_legend(ncol = 1))
# print(p_split)

## ----teardown, include = FALSE------------------------------------------------
# options(oldoption)

