#' @importFrom dplyr all_of
get_hyg <- function(x, add, id = "comid") {
  if(add && inherits(x, "sf")) {
    select(x, all_of(id))
  } else {
    NULL
  }
}

#' get node topology from edge topology (DEPRECATED)
#' @description creates a node topology table from an edge topology
#' @inheritParams get_sorted
#' @inheritParams get_tocomid
#' @param add_div data.frame containing id and toid diverted paths to add.
#' Should have id and toid fields in the first and second columns. Names
#' are not used.
#' @return data.frame containing id, fromnode, and tonode attributes or all
#' attributes provided with id, fromnode and tonode in the first three columns.
#' @export

make_node_topology <- function(x, add_div = NULL, add = TRUE) {

  warning("nhdplusTools make_node_topology is deprecated. Use hydroloom version.")

  orig_name <- names(x)[1:2]

  names(x)[1:2] <- c("id", "toid")

  x <- hydroloom::make_node_topology(x, add_div, add)

  if(add) {

    names(x)[1:2] <- orig_name[1:2]

    x

  } else {

    x <- as.data.frame(select(x, "id", "fromnode", "tonode"))

    names(x)[1] <- orig_name[1]

    x

  }
}
