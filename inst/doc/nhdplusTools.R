## ----setup, include = FALSE---------------------------------------------------
library(nhdplusTools)
local <- (Sys.getenv("BUILD_VIGNETTES") == "TRUE")
if(local) {
  cache_path <- file.path(nhdplusTools_data_dir(), "nhdpt_v_cache")
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

## ----cran, echo = TRUE, eval = FALSE------------------------------------------
#  # install.packages("nhdplusTools")

## ----install, echo = TRUE, eval = FALSE---------------------------------------
#  # install.packages("remotes")
#  # remotes::install_github("usgs-r/nhdplusTools")

## ----teardown, include=FALSE--------------------------------------------------
if(Sys.getenv("BUILD_VIGNETTES") != "TRUE") {
  unlink(work_dir, recursive = TRUE)
}


