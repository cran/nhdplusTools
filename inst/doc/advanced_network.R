## ----setup, include = FALSE---------------------------------------------------
library(nhdplusTools)
library(dplyr)

local <- (Sys.getenv("BUILD_VIGNETTES") == "TRUE")
if(local) {
  cache_path <- file.path(nhdplusTools_data_dir(), "nhdpt_v_cache")
} else {
  cache_path <- tempdir()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>",
  fig.width=4, 
  fig.height=4, 
  fig.align = "center",
  eval=local,
  cache=local,
  cache.path=(cache_path)
)

oldoption <- options(scipen = 9999,
                     "rgdal_show_exportToProj4_warnings"="none")


## ---- echo=FALSE, eval=TRUE, , fig.dim=c(6, 4)--------------------------------
print.data.frame(data.frame(ID = c(1, 2, 3), 
                            toID = c(3, 3, NA),
                            fromnode = c("N1", "N2", "N3"),
                            tonode = c("N3", "N3", "N4")), 
                 row.names = FALSE)

## ----node, fig.show="hold", out.width="45%", echo=FALSE, eval=TRUE, fig.cap="In an edge-node topology, edges are directed to nodes which are then directed to other edges. An edge-to-edge toplogy does not have intervening nodes."----
x <- c(1, 5, 3, 3)
y <- c(5, 5, 3, 1)

par(mar = c(0, 0, 0, 0))
plot(x, y, col = NA)
arrows(x[1] + 0.1, y[1] - 0.1, x[3] - 0.1, y [3] + 0.1, 0.1)
arrows(x[2] - 0.1, y[2] -0.1, x[3] + 0.1, y [3] + 0.1, 0.1)
arrows(x[3], y[3] - 0.1, x[4], y [4] + 0.1, 0.1)
text(c(2, 4, 3.15), c(4.2, 4.2, 2), c("1", "2", "3"))

par(mar = c(0, 0, 0, 0))
plot(x, y)
arrows(x[1] + 0.1, y[1] - 0.1, x[3] - 0.1, y [3] + 0.1, 0.1)
arrows(x[2] - 0.1, y[2] -0.1, x[3] + 0.1, y [3] + 0.1, 0.1)
arrows(x[3], y[3] - 0.1, x[4], y [4] + 0.1, 0.1)
text(c(2, 4, 3.1), c(4.2, 4.2, 2), c("1", "2", "3"))
text(c(1, 5, 3, 3.25), c(4.8, 4.8, 3.4, 1), c("N1", "N2", "N3", "N4"))

## ----sort, echo=FALSE, eval=TRUE, fig.cap="Smaller 'hydrosequence' values are guaranteed to be downstream of larger values along connected paths.", fig.dim=c(3, 3)----
source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))

new_hope_flowline <- get_sorted(new_hope_flowline[c("Hydroseq", "DnHydroseq")])
new_hope_flowline$sort <- seq(nrow(new_hope_flowline), 1)
par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(new_hope_flowline), col = NA)
plot(new_hope_flowline["sort"], add = TRUE, lwd = 2)

## ----lp, echo=FALSE, eval=TRUE, fig.cap="Levelpath values are constant along mainstem paths and are derived from the hydrosequence of their outlet flowline.", fig.dim=c(3, 3)----
source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))

lp <- data.frame(lp = sort(unique(new_hope_flowline$LevelPathI)))
lp$id <- seq_len(nrow(lp))
new_hope_flowline <- dplyr::left_join(new_hope_flowline, lp, by = c("LevelPathI" = "lp"))

par(mar = c(0, 0, 0, 0))
plot(sf::st_geometry(new_hope_flowline), col = NA)
plot(new_hope_flowline["id"], add = TRUE, lwd = 2)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
# Import data
source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))

# Strip the data back to the required base attributes
fpath <- get_tocomid(
  dplyr::select(new_hope_flowline, COMID, FromNode, ToNode, Divergence, FTYPE,
                AreaSqKM, LENGTHKM, GNIS_ID)
)

# Print
head(fpath <- select(sf::st_cast(fpath, "LINESTRING"), 
                     -tonode, -fromnode, -divergence, -ftype))

## -----------------------------------------------------------------------------
head(fpath <- get_sorted(fpath, split = TRUE))

## ---- echo = TRUE, fig.dim=c(3, 3)--------------------------------------------
fpath['hydrosequence'] <- seq(nrow(fpath), 1)
plot(fpath['hydrosequence'], key.pos = NULL)

## -----------------------------------------------------------------------------
# Rename and compute weight
fpath[["arbolatesum"]] <- calculate_arbolate_sum(
  dplyr::select(fpath, 
                ID = comid, toID = tocomid, length = lengthkm))

plot(sf::st_geometry(fpath), lwd = fpath$arbolatesum / 10)

## ----levelpath, fig.show="hold", out.width="45%"------------------------------
# Get levelpaths
lp <- get_levelpaths(
  dplyr::select(fpath, 
                ID = comid, toID = tocomid, 
                nameID = gnis_id, weight = arbolatesum), 
  status = FALSE, override_factor = 5)

# Print
head(fpath <- dplyr::left_join(fpath, lp, by = c("comid" = "ID")))

plot(fpath["topo_sort"], key.pos = NULL, reset = FALSE)
plot(fpath["levelpath"], key.pos = NULL)

## -----------------------------------------------------------------------------
# Invert plotting order
fpath <- dplyr::arrange(fpath, topo_sort) 

# Level Paths with more then 2 flowpaths
lp <- dplyr::group_by(fpath, levelpath) %>%
dplyr::filter(n() > 2) 

# Unique Level Path ID
lp <-  unique(lp$levelpath)

# Terminal Flowpath 
terminal_fpath <- dplyr::filter(fpath, comid %in% terminalID)

gif_file <- "levelpath.gif"

gifski::save_gif({
  for(i in 1:length(lp)) {
    lp_plot <- dplyr::filter(fpath, levelpath == lp[i])

    outlet_plot <- dplyr::filter(lp_plot, comid %in% outletID)

    plot(sf::st_geometry(fpath), lwd = 0.5, col = "grey")
    plot(sf::st_geometry(terminal_fpath), lwd = 3, col = "red", add = TRUE)
    plot(sf::st_geometry(dplyr::filter(fpath, levelpath %in% lp[1:i])), add = TRUE)
    plot(sf::st_geometry(lp_plot), col = "blue", add = TRUE)
    plot(sf::st_geometry(outlet_plot), col = "red", lwd = 1.5, add = TRUE)
  }
}, gif_file, delay = 0.5)

knitr::include_graphics(gif_file)

## -----------------------------------------------------------------------------
head(add_plus_network_attributes(dplyr::select(fpath, comid, tocomid, lengthkm, areasqkm, 
                                               nameID = gnis_id), status = TRUE))

## ----teardown, include=FALSE--------------------------------------------------
options(oldoption)

if(Sys.getenv("BUILD_VIGNETTES") != "TRUE") {
  unlink(work_dir, recursive = TRUE)
}

