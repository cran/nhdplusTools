## ----setup, include = FALSE---------------------------------------------------
library(nhdplusTools)

local <- (Sys.getenv("BUILD_VIGNETTES") == "TRUE")

if(local) {
  cache_path <- file.path(nhdplusTools_data_dir(), "data_v_cache")
} else {
  cache_path <- tempdir()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4,
  eval=local,
  cache=local,
  cache.path=(cache_path)
)

oldoption <- options(scipen = 9999)


## -----------------------------------------------------------------------------

demo_dir <- nhdplusTools_data_dir()

site_id <- "USGS-04074950"


# use dataRetrieval::get_nldi_sources() to find other nldi sources

site <- list(featureSource = "nwissite",
             featureID = "USGS-04074950")

site_feature <- get_nldi_feature(site)

upstream_network <- navigate_nldi(site, 
                                  mode = "UT", distance_km = 9999)

demo_data <- file.path(demo_dir, "data_demo.gpkg")

dataset <- subset_nhdplus(as.integer(upstream_network$UT_flowlines$nhdplus_comid), 
                          nhdplus_data = "download", # download from a service
                          output_file = demo_data, # write the data to disk
                          return_data = TRUE, # return the data rather
                          flowline_only = FALSE, overwrite = TRUE)

names(dataset)

sapply(dataset, nrow)

old_par <- par(mar = c(0, 0, 0, 0))

plot_nhdplus(outlets = list(featureSource = "nwissite",
                            featureID = "USGS-04074950"), 
             nhdplus_data = demo_data, flowline_only = TRUE)

## -----------------------------------------------------------------------------

basin <- get_nldi_basin(site)

subset <- get_nhdplus(AOI = basin, realization = "flowline")

par(mar = c(0, 0, 0, 0))

plot(sf::st_geometry(basin))
plot(sf::st_geometry(subset), col = "blue", add = TRUE)


## -----------------------------------------------------------------------------

wolf_huc <- get_huc(basin, type = 'huc04')

nrow(wolf_huc)

# it straddles hucs? Not really.

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(basin), col = "grey")
plot(sf::st_geometry(wolf_huc), add = TRUE)

wolf_huc <- get_huc(site_feature, type = "huc04")

nrow(wolf_huc)

# better!!
par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(wolf_huc))
plot(sf::st_geometry(basin), col = "grey", add = TRUE)


## -----------------------------------------------------------------------------

outdir <- file.path(nhdplusTools_data_dir(), "hr_access_demo")

dir.create(outdir)

download_dir <- download_nhdplushr(outdir, wolf_huc$huc4)

list.files(download_dir)


## -----------------------------------------------------------------------------

nhdplushr <- get_nhdplushr(
  download_dir, 
  layers = c("NHDFlowline", "NHDPlusCatchment", "NHDWaterbody", 
             "NHDArea", "NHDLine", "NHDPlusSink", "NHDPlusWall", 
             "NHDPoint", "NHDPlusBurnWaterbody", "NHDPlusBurnLineEvent",
             "HYDRO_NET_Junctions", "WBDHU2", "WBDHU4","WBDHU6", 
             "WBDHU8", "WBDHU10", "WBDHU12", "WBDLine"), 
  check_terminals = TRUE)

sapply(nhdplushr, nrow)


## -----------------------------------------------------------------------------

nhdplushr <- get_hr_data(list.files(download_dir, pattern = ".gdb", full.names = TRUE),
                         layer = "NHDFlowline", rename = TRUE)

names(nhdplushr)

# Great Lakes coast are flowlines -- remove for visuals
gl_coast <- c(get_DM(nhdplushr, 60002700000331),
              get_DM(nhdplushr, 60002700049157))

plot_data <- dplyr::filter(nhdplushr, FCODE != 56600 & StreamOrde > 2 & !COMID %in% gl_coast)

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(wolf_huc))
plot(sf::st_geometry(basin), col = "grey", add = TRUE)
plot(sf::st_geometry(plot_data), lwd = plot_data$StreamOrde / 3, col = "blue", add = TRUE)

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(basin), col = "grey")
plot(sf::st_geometry(plot_data), lwd = plot_data$StreamOrde / 2, col = "blue", add = TRUE)

## -----------------------------------------------------------------------------

potential_matches <- get_flowline_index(nhdplushr, 
                                        points = site_feature, 
                                        search_radius = units::as_units(500, "m"),
                                        max_matches = 5)

potential_matches

site_meta <- dataRetrieval::readNWISsite(gsub("USGS-", "", site_feature$identifier))

sqmi_to_sqkm <- 2.59

da_df <- data.frame(id = 1, drainagearea = site_meta$drain_area_va * sqmi_to_sqkm)

# uses nearest drainage area match
disambiguate_flowline_indexes(potential_matches,
                              flowpath = dplyr::select(nhdplushr, COMID, TotDASqKM),
                              hydro_location = da_df)


## -----------------------------------------------------------------------------

outdir <- file.path(nhdplusTools_data_dir(), "nhd_access_demo")

dir.create(outdir)

download_dir <- download_nhd(outdir, wolf_huc$huc4)

list.files(download_dir)

nhd_gdb <- list.files(download_dir, pattern = ".gdb", full.names = TRUE)

sf::st_layers(nhd_gdb)

nhd_fline <- sf::read_sf(nhd_gdb, "NHDFlowline")


## -----------------------------------------------------------------------------

sub_3dhp <- get_3dhp(basin, type = "flowline")

plot(sub_3dhp)

sub_3dhp <- st_compatibalize(sub_3dhp, nhd_fline)

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(basin), col = "grey")
plot(sf::st_geometry(sf::st_zm(sub_3dhp)), col = "skyblue", lwd = 2.5, add = TRUE)
plot(sf::st_geometry(sf::st_zm(nhd_fline)), col = "blue", lwd = 0.1, add = TRUE)

## -----------------------------------------------------------------------------

unique(discover_geoconnex_reference()[c("id", "title")])


## -----------------------------------------------------------------------------

wolf_gages <- get_geoconnex_reference(basin, type = "gages")

geoconnex_gage <- dplyr::filter(wolf_gages, 
                                provider_id == site_feature$identifier)

wolf_mainstems <- get_geoconnex_reference(basin, type = "mainstems")

wolf_mainstem <- dplyr::filter(wolf_mainstems, uri == geoconnex_gage$mainstem_uri)

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(basin), col = "grey")
plot(sf::st_geometry(wolf_mainstems), lwd = 0.5, col = "blue", add = TRUE)
plot(sf::st_geometry(wolf_mainstem), lwd = 3, col = "blue", add = TRUE)
plot(sf::st_geometry(wolf_gages), pch = 2, cex = 0.75, add = TRUE)
plot(sf::st_geometry(geoconnex_gage), pch = 17, cex = 1.5, add = TRUE)

## -----------------------------------------------------------------------------

wbd_dir <- file.path(nhdplusTools_data_dir(), "wbd_access_demo")
  
wbd_out <- download_wbd(wbd_dir)

# zip::unzip doesn't always work
if(length(wbd_out == 0)) {
  f <- list.files(wbd_dir, pattern = ".zip", full.names = TRUE)  
  
  utils::unzip(f, exdir = wbd_dir)
}

wbd_gdb <- list.files(wbd_dir, pattern = ".gdb", full.names = TRUE)

sf::st_layers(wbd_gdb)


## ----teardown, include=FALSE--------------------------------------------------
options(oldoption)
par(old_par)
if(Sys.getenv("BUILD_VIGNETTES") != "TRUE") {
  unlink(work_dir, recursive = TRUE)
}

