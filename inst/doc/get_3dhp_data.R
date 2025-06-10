## ----setup, include = FALSE---------------------------------------------------
library(nhdplusTools)

local <- (Sys.getenv("BUILD_VIGNETTES") == "TRUE")

if(local) {
  cache_path <- file.path(nhdplusTools_data_dir(), "3dhp_v_cache")
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


## ----message=FALSE------------------------------------------------------------

  comid <- "937070225"

  point <- c(-89.441, 43.487) |>
    sf::st_point() |>
    sf::st_sfc(crs = 4326) |> sf::st_sf()

  dm <- c('https://geoconnex.us/ref/mainstems/323742',
          'https://geoconnex.us/ref/mainstems/312091')

  flowline <- get_3dhp(point, type = "flowline", buffer = 10)
  
  pc <- function(x) sf::st_geometry(sf::st_transform(x, 3857))

  plot_nhdplus(bbox = sf::st_bbox(flowline), plot_config = list(flowline = list(col = NULL)), zoom = 14)
  plot(pc(flowline), add = TRUE, col = "darkblue", lwd = 2)
  plot(pc(point), pch = "O", add = TRUE)

## ----message=FALSE------------------------------------------------------------
  basin <- dataRetrieval::findNLDI(comid = comid, find = "basin")

  network <- get_3dhp(basin$basin, type = "flowline")
  water <- get_3dhp(basin$basin, type = "waterbody")
  hydrolocation <- get_3dhp(basin$basin, type = "hydrolocation - reach code, external connection")
  
  down_mains <- get_3dhp(ids = dm, type = "flowline")
  
  old_par <- par(mar = c(0, 0, 0, 0))
  plot_nhdplus(bbox = sf::st_bbox(basin$basin), flowline_only = TRUE,
               plot_config = list(flowline = list(col = NULL)), zoom = 10)
  plot(pc(basin$basin), lwd = 2, add = TRUE)
  plot(pc(network), lwd = 0.5, add = TRUE)
  plot(pc(water), lwd = 0.5, border = "skyblue", col = "lightblue", add = TRUE)
  plot(pc(hydrolocation), pch = "o", col = "#80808026", add = TRUE)
  plot(pc(down_mains), lwd = 3, col = "blue", add = TRUE)
  par(old_par)

  old_par <- par(mar = c(0, 0, 0, 0))
  plot(pc(down_mains))
  plot(pc(basin$basin), add = TRUE)
  par(old_par)

## -----------------------------------------------------------------------------

  reachcode <- "07070004002889"

  hydrolocation <- get_3dhp(universalreferenceid = reachcode,
                            type = "hydrolocation - reach code, external connection")

  hydrolocation

  mainstem_points <- get_3dhp(ids = hydrolocation$mainstemid, type = "hydrolocation - reach code, external connection")
  mainstem_lines <- get_3dhp(ids = hydrolocation$mainstemid, type = "flowline")

  old_par <- par(mar = c(0, 0, 0, 0))
  plot(pc(hydrolocation), pch = "o", col = "#808080BF")
  plot(pc(mainstem_lines), lwd = 2, col = "blue", add = TRUE)
  par(old_par)
  
  old_par <- par(mar = c(0, 0, 0, 0))
  plot_nhdplus(bbox = sf::st_bbox(basin$basin), flowline_only = TRUE,
               plot_config = list(flowline = list(col = NULL)), zoom = 10)
  plot(pc(mainstem_lines), lwd = 2, col = "blue", add = TRUE)
  plot(pc(mainstem_points), pch = "o", col = "#808080BF", add = TRUE)
  par(old_par)

## ----teardown, include=FALSE--------------------------------------------------
options(oldoption)
par(old_par)
if(Sys.getenv("BUILD_VIGNETTES") != "TRUE") {
  unlink(work_dir, recursive = TRUE)
}

