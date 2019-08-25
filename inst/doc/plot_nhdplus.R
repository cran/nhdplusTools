## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4
)
oldoption <- options(scipen = 9999)
options(scipen = 9999)

## ----data_dir_setup, echo=FALSE, include=FALSE---------------------------
temp_dir <- tempdir()
dir.create(temp_dir)

## ----get data------------------------------------------------------------
library(sf)
library(nhdplusTools)
nwissite <- list(featureSource = "nwissite", 
                     featureID = "USGS-05428500")

flowline <- navigate_nldi(nwissite, 
                          mode = "upstreamTributaries", 
                          data_source = "")

nhdplus <- subset_nhdplus(comids = flowline$nhdplus_comid,
                          output_file = file.path(temp_dir, "nhdplus.gpkg"),
                          nhdplus_data = "download",
                          overwrite = TRUE)

flowline <- read_sf(nhdplus, "NHDFlowline_Network")

upstream_nwis <- navigate_nldi(nwissite,
                               mode = "upstreamTributaries",
                               data_source = "nwissite")

basin <- get_nldi_basin(nwissite)

## ----introspect----------------------------------------------------------
st_layers(nhdplus)
names(flowline)
names(upstream_nwis)
names(basin)
class(st_geometry(flowline))
class(st_geometry(upstream_nwis))
class(st_geometry(basin))

## ----bbox, message=FALSE-------------------------------------------------
library(sp)

sf_bbox <- st_bbox(basin)
sf_bbox
class(sf_bbox)

sp_bbox <- sp::bbox(sf::as_Spatial(basin))
sp_bbox
class(sp_bbox)

# Or without the sp::bbox
sp_bbox <- matrix(sf_bbox, 
                  byrow = FALSE, 
                  ncol = 2, 
                  dimnames = list(c("x", "y"), 
                                  c("min", "max")))
sp_bbox

ggmap_bbox <- setNames(sf_bbox, c("left", "bottom", "right", "top"))
ggmap_bbox

## ----plot----------------------------------------------------------------
prep_layer <- function(x) st_geometry(st_transform(x, 3857))

prettymapr::prettymap({
  rosm::osm.plot(sp_bbox, type = "stamenwatercolor", quiet = TRUE)
  
  plot(prep_layer(basin), 
       lwd = 2, add = TRUE)
  
  plot(prep_layer(flowline), 
       lwd = 1.5, col = "deepskyblue", add = TRUE)
  
  plot(prep_layer(dplyr::filter(flowline, streamorde > 2)), 
       lwd = 3, col = "darkblue", add = TRUE)
  
  us_nwis_layer <- prep_layer(upstream_nwis)
  
  plot(us_nwis_layer, 
       pch = 17, cex = 1.5, col = "yellow", add = TRUE)
  
  label_pos <- st_coordinates(us_nwis_layer)
  
  text(label_pos[,1],label_pos[,2], 
       upstream_nwis$identifier, 
       adj = c(-0.2, 0.5), cex = 0.7)
  
}, drawarrow = TRUE)

## ----ggmap---------------------------------------------------------------
library(ggmap)
library(ggplot2)

upstream_nwis[c("lon", "lat")] <- sf::st_coordinates(upstream_nwis)

basemap_toner <- get_map(source = "stamen", maptype = "toner", location = ggmap_bbox, zoom = 11)
basemap_terrain <- get_map(source = "stamen", maptype = "terrain-lines", location = ggmap_bbox, zoom = 11)
toner_map <- ggmap(basemap_toner)
terrain_map <- ggmap(basemap_terrain)

toner_map

terrain_map + geom_sf(data = basin,
                        inherit.aes = FALSE,
                        color = "black", fill = NA) + 
  geom_sf(data = flowline,
          inherit.aes = FALSE,
          color = "deepskyblue") +
  geom_sf(data = dplyr::filter(flowline, streamorde > 2),
          inherit.aes = FALSE,
          color = "darkblue") +
  geom_sf(data = upstream_nwis, inherit.aes = FALSE, color = "red") + 
  geom_text(data = upstream_nwis, aes(label = identifier, x = lon, y = lat),
            hjust = 0, size=2.5, nudge_x = 0.02, col = "black")

## ----teardown, include=FALSE---------------------------------------------
options(oldoption)

