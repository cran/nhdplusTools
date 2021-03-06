# Primary HydroShare Data Resource
vaa_hydroshare <-
  'https://www.hydroshare.org/resource/6092c8a62fac45be97a09bfd0b0bf726/data/contents/nhdplusVAA.fst'

nhdplusTools_env <- new.env()

# NHDPlus Attributes
COMID <- "COMID"
FEATUREID <- "FEATUREID"
Hydroseq <- "Hydroseq"
DnHydroseq <- "DnHydroseq"
DnMinorHyd <- "DnMinorHyd"
LevelPathI <- "LevelPathI"
DnLevelPat <- "DnLevelPat"
ToNode <- "ToNode"
FromNode <- "FromNode"
TotDASqKM <- "TotDASqKM"
AreaSqKM <- "AreaSqKM"
LENGTHKM <- "LENGTHKM"
Pathlength <- "Pathlength"
StreamCalc <- "StreamCalc"
StreamOrde <- "StreamOrde"
TerminalFl <- "TerminalFl"
Divergence <- "Divergence"
TerminalPa <- "TerminalPa"
StartFlag <- "StartFlag"
FTYPE <- "FTYPE"
FCODE <- "FCODE"
FromMeas <- "FromMeas"
ToMeas <- "ToMeas"
REACHCODE <- "REACHCODE"
REACH_meas <- "REACH_meas"
HUC12 <- "HUC12"
TOHUC <- "TOHUC"
ReachCode <- "ReachCode"
VPUID <- "VPUID"
toCOMID = "toCOMID"


# List of input names that should be changed to replacement names
nhdplus_attributes <- list(
  COMID = COMID, NHDPlusID = COMID,
  FEATUREID = FEATUREID,
  Hydroseq = Hydroseq, HydroSeq = Hydroseq,
  DnHydroseq = DnHydroseq, DnHydroSeq = DnHydroseq,
  DnMinorHyd = DnMinorHyd,
  LevelPathI = LevelPathI,
  DnLevelPat = DnLevelPat,
  ToNode = ToNode,
  FromNode = FromNode,
  TotDASqKM = TotDASqKM, TotDASqKm = TotDASqKM,
  AreaSqKM = AreaSqKM, AreaSqKm = AreaSqKM,
  LENGTHKM = LENGTHKM, LengthKM = LENGTHKM,
  Pathlength = Pathlength, PathLength = Pathlength,
  StreamCalc = StreamCalc,
  StreamOrde = StreamOrde,
  TerminalFl = TerminalFl,
  Divergence = Divergence,
  TerminalPa = TerminalPa,
  StartFlag = StartFlag,
  FTYPE = FTYPE, FType = FTYPE,
  FCODE = FCODE, FCode = FCODE,
  FromMeas = FromMeas,
  ToMeas = ToMeas,
  REACHCODE = REACHCODE, ReachCode = REACHCODE,
  REACH_meas = REACH_meas,
  HUC12 = HUC12,
  TOHUC = TOHUC,
  toCOMID = toCOMID)

.data <- . <- NULL

assign("nhdplus_attributes", nhdplus_attributes, envir = nhdplusTools_env)

assign("geoserver_root", "https://labs.waterdata.usgs.gov/geoserver/",
       envir = nhdplusTools_env)

assign("prepare_nhdplus_attributes",
       c("COMID", "LENGTHKM", "TerminalFl",
         "FromNode", "ToNode", "TotDASqKM",
         "StartFlag", "StreamOrde", "StreamCalc",
         "TerminalPa", "Pathlength", "Divergence", "Hydroseq",
         "LevelPathI"),
       envir = nhdplusTools_env)

assign("split_flowlines_attributes",
       c("COMID", "toCOMID", "LENGTHKM"),
       envir = nhdplusTools_env)

assign("collapse_flowlines_attributes",
       c("COMID", "toCOMID", "LENGTHKM", "LevelPathI", "Hydroseq"),
       envir = nhdplusTools_env)

assign("reconcile_collapsed_flowlines_attributes",
       c("COMID", "toCOMID", "LENGTHKM", "LevelPathI", "Hydroseq"),
       envir = nhdplusTools_env)

assign("get_UT_attributes",
       c("COMID", "Pathlength", "LENGTHKM", "Hydroseq",
         "LevelPathI", "DnHydroseq"),
       envir = nhdplusTools_env)

assign("get_UM_attributes",
       c("COMID", "Pathlength", "LevelPathI",
         "Hydroseq"),
       envir = nhdplusTools_env)

assign("get_DM_attributes",
       c("COMID", "Pathlength", "LENGTHKM",
         "LevelPathI", "DnLevelPat",
         "DnHydroseq", "Hydroseq"),
       envir = nhdplusTools_env)

assign("get_DM_nolength_attributes",
       c("COMID",
         "LevelPathI", "DnLevelPat",
         "DnHydroseq", "Hydroseq"),
       envir = nhdplusTools_env)

assign("get_DD_attributes",
       c("COMID", "Pathlength", "LENGTHKM",
         "LevelPathI", "DnLevelPat",
         "DnHydroseq", "Hydroseq", "DnMinorHyd"),
       envir = nhdplusTools_env)

assign("get_flowline_index_attributes",
       c("COMID", "REACHCODE", "ToMeas", "FromMeas"),
       envir = nhdplusTools_env)

assign("get_levelpaths_attributes",
       c("ID", "toID", "nameID", "weight"),
       envir = nhdplusTools_env)

assign("get_streamorder_attributes",
       c("ID", "toID"),
       envir = nhdplusTools_env)

assign("get_pfaf_attributes",
       c("ID", "toID", "totda", "outletID", "topo_sort", "levelpath"),
       envir = nhdplusTools_env)

assign("make_standalone_tonode_attributes",
       c("COMID", "ToNode", "FromNode", "TerminalFl", "Hydroseq", "TerminalPa",
         "LevelPathI", "FTYPE"), envir = nhdplusTools_env)

assign("make_standalone_tocomid_attributes",
       c("COMID", "toCOMID", "Hydroseq", "TerminalPa",
         "LevelPathI", "FTYPE"), envir = nhdplusTools_env)

assign("get_waterbody_index_waterbodies_attributes",
       c("COMID"), envir = nhdplusTools_env)

assign("get_waterbody_index_flines_attributes",
       c("COMID", "WBAREACOMI", "Hydroseq"),
       envir = nhdplusTools_env)

assign("disambiguate_flowline_indexes_attributes",
       c("id", "COMID", "REACHCODE", "REACH_meas", "offset"),
       envir = nhdplusTools_env)

# assigned here for record keeping. Used as a status counter in apply functions.
assign("cur_count", 0, envir = nhdplusTools_env)

check_names <- function(x, function_name) {
  x <- align_nhdplus_names(x)
  names_x <- names(x)
  expect_names <- get(paste0(function_name, "_attributes"),
                      envir = nhdplusTools_env)
  if ( !all(expect_names %in% names_x)) {
    stop(paste0("Missing some required attributes in call to: ",
                function_name, ". Expected: ",
                paste(expect_names[which(!(expect_names %in%
                                             names_x))],
                      collapse = ", "), " or NHDPlusHR equivalents."))
  }
  return(x)
}

default_nhdplus_path <- "../NHDPlusV21_National_Seamless.gdb"

assign("default_nhdplus_path", default_nhdplus_path, envir = nhdplusTools_env)

nhdhr_bucket <- "https://prd-tnm.s3.amazonaws.com/"
nhdhr_file_list <- "?prefix=StagedProducts/Hydrography/NHDPlusHR/Beta/GDB/"

assign("nhdhr_bucket", nhdhr_bucket, envir = nhdplusTools_env)
assign("nhdhr_file_list", nhdhr_file_list, envir = nhdplusTools_env)

assign("nhdpt_dat_dir",
       tools::R_user_dir("usgs_r/nhdplusTools"),
       envir = nhdplusTools_env)

#' get or set nhdplusTools data directory
#' @description if left unset, will return the user data dir
#' as returned by \link[tools]{R_user_dir} for this package.
#' @param dir path of desired data directory
#' @return character path of data directory (silent when setting)
#' @importFrom tools R_user_dir
#' @export
#' @examples
#' nhdplusTools_data_dir()
#'
#' nhdplusTools_data_dir("demo")
#'
#' nhdplusTools_data_dir(tools::R_user_dir("usgs_r/nhdplusTools"))
#'
nhdplusTools_data_dir <- function(dir = NULL) {

  if(is.null(dir)) {
    return(get("nhdpt_dat_dir", envir = nhdplusTools_env))
  } else {
    assign("nhdpt_dat_dir",
           dir,
           envir = nhdplusTools_env)
    return(invisible(get("nhdpt_dat_dir", envir = nhdplusTools_env)))
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste(strwrap(
    "USGS Support Package:
    https://owi.usgs.gov/R/packages.html#support"),
    collapse = "\n"))
  nhdplus_path(default_nhdplus_path, warn = FALSE)
}

#' @title NHDPlus Data Path
#' @description Allows specification of a custom path to a source dataset.
#' Typically this will be the national seamless dataset in
#' geodatabase or geopackage format.
#' @param path character path ending in .gdb or .gpkg
#' @param warn boolean controls whether warning an status messages are printed
#' @return 0 (invisibly) if set successfully, character path if no input.
#' @export
#' @examples
#' nhdplus_path("/data/NHDPlusV21_National_Seamless.gdb")
#'
#' nhdplus_path("/data/NHDPlusV21_National_Seamless.gdb", warn=FALSE)
#'
#' nhdplus_path()
#'
nhdplus_path <- function(path = NULL, warn = FALSE) {
  if (!is.null(path)) {

    assign("nhdplus_data", path, envir = nhdplusTools_env)

    if (warn) {
      warning("Path does not exist.")
    }

    if (nhdplus_path() == path) {
      invisible(0)
    }
  } else {
      return(get("nhdplus_data", envir = nhdplusTools_env))
  }
}


#' @title Align NHD Dataset Names
#' @description this function takes any NHDPlus dataset and aligns the attribute names with those used in nhdplusTools.
#' @param x a \code{sf} object of nhdplus flowlines
#' @return data.frame renamed \code{sf} object
#' @export
#' @examples
#' source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))
#'
#' names(new_hope_flowline)
#'
#' names(new_hope_flowline) <- tolower(names(new_hope_flowline))
#'
#' new_hope_flowline <- align_nhdplus_names(new_hope_flowline)
#'
#' names(new_hope_flowline)
#'
align_nhdplus_names <- function(x){

  attribute_names <- get("nhdplus_attributes", envir = nhdplusTools_env)

  # get into correct case
  good_names <- unique(unlist(do.call(rbind, attribute_names))[,1])

  new_names <- old_names <- names(x)

  matched <- match(toupper(names(x)), toupper(good_names))
  replacement_names <- as.character(good_names[matched[which(!is.na(matched))]])

  new_names[which(toupper(old_names) %in% toupper(good_names))] <- replacement_names
  names(x) <- new_names

  # rename to match package
  new_names <- old_names <- names(x)

  matched <- match(names(x), names(attribute_names))
  replacement_names <- as.character(attribute_names[matched[which(!is.na(matched))]])

  new_names[which(old_names %in% names(attribute_names))] <- replacement_names

  names(x) <- new_names

  if("GridCode" %in% names(x) & !"FeatureID" %in% names(x) & !"FEATUREID" %in% names(x))
    names(x)[which(names(x) == "COMID")] <- "FEATUREID"

  return(x)

}


drop_geometry <- function(x) {
  if("sf" %in% class(x)) {
    sf::st_drop_geometry(x)
  } else {
    x
  }
}

get_cl <- function(cl) {
  if(!is.null(cl)) {
    if(!requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel required if using cores input")
    }
    if(is.numeric(cl)) {
      cl <- parallel::makeCluster(cl)
    } else {
      if(!"cluster" %in% class(cl)) {
        stop("cores must be numeric or a cluster object")
      }
    }
  }
  return(cl)
}
