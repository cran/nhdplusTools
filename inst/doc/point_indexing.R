## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4,
  eval=nzchar(Sys.getenv("BUILD_VIGNETTES"))
)
oldoption <- options(scipen = 9999)
options(scipen = 9999)

## ----nhdplus_path_setup, echo=FALSE, include=FALSE----------------------------
#  library(dplyr, warn.conflicts = FALSE)
#  
#  temp_dir <- tempdir()
#  dir.create(temp_dir)
#  
#  file.copy(system.file("extdata/sample_natseamless.gpkg",
#                        package = "nhdplusTools"),
#            file.path(temp_dir, "natseamless.gpkg"))

## ----nhdplus_path, echo=TRUE--------------------------------------------------
#  library(nhdplusTools)
#  
#  nhdplus_path(file.path(temp_dir, "natseamless.gpkg"))
#  
#  flowlines <- sf::read_sf(nhdplus_path(), "NHDFlowline_Network")
#  gages <- sf::read_sf(nhdplus_path(), "Gage")

## ----get_indexes--------------------------------------------------------------
#  geom_col <- attr(gages, "sf_column")
#  indexes <- get_flowline_index(flowlines, gages[[geom_col]], search_radius = 0.1)

## ----analyze_index------------------------------------------------------------
#  p_match <- 100 * length(which(indexes$COMID %in% gages$FLComID)) / nrow(gages)
#  paste0(round(p_match, digits = 1),
#         "% were found to match the COMID in the NHDPlus gages layer")
#  
#  p_match <- 100 * length(which(indexes$REACHCODE %in% gages$REACHCODE)) / nrow(gages)
#  paste0(round(p_match, digits = 1),
#         "% were found to match the REACHCODE in the NHDPlus gages layer")
#  
#  matched <- cbind(indexes,
#                   dplyr::select(sf::st_set_geometry(gages, NULL),
#                                 REACHCODE_ref = REACHCODE,
#                                 COMID_ref = FLComID,
#                                 REACH_meas_ref = Measure)) %>%
#    dplyr::filter(REACHCODE == REACHCODE_ref) %>%
#    dplyr::mutate(REACH_meas_diff = REACH_meas - REACH_meas_ref)
#  
#  hist(matched$REACH_meas_diff, breaks = 100,
#       main = "Difference in measure for gages matched to the same reach.")
#  
#  round(quantile(matched$REACH_meas_diff,
#                 probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)),
#        digits = 2)

## ----get_indexes_precise------------------------------------------------------
#  indexes <- get_flowline_index(flowlines,
#                                gages[[geom_col]],
#                                search_radius = 0.1,
#                                precision = 10)

## ----analyze_inde_precise-----------------------------------------------------
#  p_match <- 100 * length(which(indexes$COMID %in% gages$FLComID)) / nrow(gages)
#  paste0(round(p_match, digits = 1),
#         "% were found to match the COMID in the NHDPlus gages layer")
#  
#  p_match <- 100 * length(which(indexes$REACHCODE %in% gages$REACHCODE)) / nrow(gages)
#  paste0(round(p_match, digits = 1),
#         "% were found to match the REACHCODE in the NHDPlus gages layer")
#  
#  matched <- cbind(indexes,
#                   dplyr::select(sf::st_set_geometry(gages, NULL),
#                                 REACHCODE_ref = REACHCODE,
#                                 COMID_ref = FLComID,
#                                 REACH_meas_ref = Measure)) %>%
#    dplyr::filter(REACHCODE == REACHCODE_ref) %>%
#    dplyr::mutate(REACH_meas_diff = REACH_meas - REACH_meas_ref)
#  
#  hist(matched$REACH_meas_diff, breaks = 100,
#       main = "Difference in measure for gages matched to the same reach.")
#  
#  round(quantile(matched$REACH_meas_diff,
#                 probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), digits = 2)

## ----teardown, include=FALSE--------------------------------------------------
#  options(oldoption)

