## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4,
  eval=nzchar(Sys.getenv("BUILD_VIGNETTES"))
)
oldoption <- options(scipen = 9999)
options(scipen = 9999)

## ----tldr----------------------------------------------------------------
#  # Uncomment to install!
#  # install.packages("nhdplusTools")
#  
#  library(nhdplusTools)
#  library(sf)
#  
#  start_point <- st_sfc(st_point(c(-89.362239, 43.090266)), crs = 4269)
#  start_comid <- discover_nhdplus_id(start_point)
#  
#  flowline <- navigate_nldi(list(featureSource = "comid",
#                                 featureID = start_comid),
#                            mode = "upstreamTributaries",
#                            data_source = "")
#  
#  subset_gpkg <-subset_nhdplus(comids = flowline$nhdplus_comid,
#                               output_file = tempfile(fileext = ".gpkg"),
#                               nhdplus_data = "download")
#  
#  flowline <- sf::read_sf(subset_gpkg, "NHDFlowline_Network")
#  catchment <- sf::read_sf(subset_gpkg, "CatchmentSP")
#  waterbody <- sf::read_sf(subset_gpkg, "NHDWaterbody")
#  
#  plot(sf::st_geometry(flowline), col = "blue")
#  plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
#  plot(sf::st_geometry(catchment), add = TRUE)
#  plot(sf::st_geometry(waterbody), col = rgb(0, 0, 1, alpha = 0.5), add = TRUE)

## ----install, echo = TRUE, eval = FALSE----------------------------------
#  # install.packages("devtools")
#  # devtools::install_github("usgs-r/nhdplusTools")

## ----load----------------------------------------------------------------
#  library(nhdplusTools)

## ----nhdplus_path_setup, echo=FALSE, include=FALSE-----------------------
#  temp_dir <- tempdir()
#  
#  dir.create(temp_dir)
#  
#  file.copy(system.file("extdata/sample_natseamless.gpkg",
#                        package = "nhdplusTools"),
#            file.path(temp_dir, "natseamless.gpkg"))

## ----nhdplus_path, echo=TRUE---------------------------------------------
#  nhdplus_path(file.path(temp_dir, "natseamless.gpkg"))
#  
#  nhdplus_path()

## ----stage_national_data-------------------------------------------------
#  staged_data <- stage_national_data(output_path = tempdir())
#  
#  str(staged_data)

## ----staged_data---------------------------------------------------------
#  flowline <- readRDS(staged_data$flowline)
#  names(flowline)[1:10]
#  
#  library(sf)
#  plot(sf::st_geometry(flowline))

## ----nhdplushr_secret, echo=FALSE, include=FALSE-------------------------
#    nhd_dir <- tempdir()
#    temp_file <- file.path(nhd_dir, "temp.zip")
#    download.file("https://usgs-r.github.io/nhdplusTools/data/03_sub.zip",
#              temp_file)
#    unzip(temp_file, exdir = nhd_dir)

## ----nhdplus_hr----------------------------------------------------------
#  download_nhdplushr(nhd_dir = "download_dir",
#                     hu_list = c("0101"), # can mix hu02 and hu04 codes.
#                     download_files = FALSE) # TRUE will download files.
#  
#  hr_data <- get_nhdplushr(nhd_dir,
#                           out_gpkg = file.path(nhd_dir, "nhd_hr.gpkg"))
#  (layers <- st_layers(hr_data))
#  unlink(hr_data)
#  
#  hr_data <- get_nhdplushr(nhd_dir,
#                           out_gpkg = file.path(nhd_dir, "nhd_hr.gpkg"),
#                           layers = NULL)
#  (layers <- st_layers(hr_data))

## ----point---------------------------------------------------------------
#  lon <- -89.362239
#  lat <- 43.090266
#  
#  start_point <- sf::st_sfc(sf::st_point(c(lon, lat)),
#                            crs = 4269)
#  
#  plot(sf::st_geometry(flowline))
#  plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)

## ----discover_nhdplus_id-------------------------------------------------
#  start_comid <- discover_nhdplus_id(start_point)
#  start_comid

## ----discover_nldi-------------------------------------------------------
#  discover_nldi_sources()$source
#  
#  nldi_feature <- list(featureSource = "comid", featureID = start_comid)
#  discover_nldi_navigation(nldi_feature)

## ----navigate_nldi-------------------------------------------------------
#  
#  flowline_nldi <- navigate_nldi(nldi_feature,
#                                 mode = "upstreamTributaries",
#                                 data_source = "")
#  
#  plot(sf::st_geometry(flowline), lwd = 3, col = "black")
#  plot(sf::st_geometry(flowline_nldi), lwd = 1, col = "red", add = TRUE)

## ----subset_nhdplus_download---------------------------------------------
#  output_file_download <- file.path(temp_dir, "subset_download.gpkg")
#  
#  output_file_download <-subset_nhdplus(comids = flowline_nldi$nhdplus_comid,
#                               output_file = output_file_download,
#                               nhdplus_data = "download")
#  
#  sf::st_layers(output_file_download)
#  
#  flowline_download <- sf::read_sf(file.path(temp_dir, "subset_download.gpkg"),
#                                   "NHDFlowline_Network")
#  
#  plot(sf::st_geometry(dplyr::filter(flowline_download,
#                                     streamorde > 2)),
#       lwd = 7, col = "darkgrey")
#  plot(sf::st_geometry(flowline_nldi),
#       lwd = 3, col = "red", add = TRUE)

## ----nldi_nwissite-------------------------------------------------------
#  nldi_feature <- list(featureSource = "nwissite", featureID = "USGS-05428500")
#  
#  flowline_nldi <- navigate_nldi(nldi_feature,
#                                 mode = "upstreamTributaries",
#                                 data_source = "")
#  
#  output_file_nwis <- file.path(temp_dir, "subset_download_nwis.gpkg")
#  
#  output_file_nwis <-subset_nhdplus(comids = flowline_nldi$nhdplus_comid,
#                                    output_file = output_file_nwis,
#                                    nhdplus_data = "download")
#  
#  sf::st_layers(output_file_download)
#  
#  flowline_nwis <- sf::read_sf(output_file_nwis,
#                                   "NHDFlowline_Network")
#  
#  upstream_nwis <- navigate_nldi(nldi_feature,
#                                 mode = "upstreamTributaries",
#                                 data_source = "nwissite")
#  
#  plot(sf::st_geometry(flowline_nwis),
#       lwd = 3, col = "blue")
#  plot(sf::st_geometry(upstream_nwis),
#       cex = 1, lwd = 2, col = "red", add = TRUE)

## ----get_UT--------------------------------------------------------------
#  UT_comids <- get_UT(flowline, start_comid)
#  UT_comids

## ----plot_fline_subset---------------------------------------------------
#  plot(sf::st_geometry(flowline))
#  plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
#  plot(sf::st_geometry(dplyr::filter(flowline, COMID %in% UT_comids)),
#       add=TRUE, col = "red", lwd = 2)

## ----subset_nhdplus------------------------------------------------------
#  output_file <- file.path(temp_dir, "subset.gpkg")
#  
#  output_file <-subset_nhdplus(comids = UT_comids,
#                               output_file = output_file,
#                               nhdplus_data = nhdplus_path())
#  
#  sf::st_layers(output_file)

## ----plot_result---------------------------------------------------------
#  catchment <- sf::read_sf(output_file, "CatchmentSP")
#  waterbody <- sf::read_sf(output_file, "NHDWaterbody")
#  
#  plot(sf::st_geometry(flowline))
#  plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
#  plot(sf::st_geometry(dplyr::filter(flowline, COMID %in% UT_comids)),
#       add=TRUE, col = "red", lwd = 2)
#  plot(sf::st_geometry(catchment), add = TRUE)
#  plot(sf::st_geometry(waterbody), col = rgb(0, 0, 1, alpha = 0.5), add = TRUE)

## ----indexing------------------------------------------------------------
#  get_flowline_index(flowline, start_point)

## ----index_list----------------------------------------------------------
#  gage <- sf::read_sf(output_file, "Gage")
#  
#  get_flowline_index(flowline, sf::st_geometry(gage), precision = 10)

## ----teardown, include=FALSE---------------------------------------------
#  options(oldoption)

