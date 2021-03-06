context("calculate network attributes")

test_that("total drainage area works", {
  source(system.file("extdata", "walker_data.R", package = "nhdplusTools"))

  catchment_area <- prepare_nhdplus(walker_flowline, 0, 0,
                                    purge_non_dendritic = FALSE, warn = FALSE)

  catchment_area <- select(walker_flowline, COMID, AreaSqKM) %>%
    left_join(catchment_area, by = "COMID") %>%
    select(ID = COMID, toID = toCOMID, area = AreaSqKM)

  new_da <- calculate_total_drainage_area(catchment_area)

  catchment_area$totda <- new_da
  catchment_area$nhdptotda <- walker_flowline$TotDASqKM

  expect(mean(abs(catchment_area$totda - catchment_area$nhdptotda)) < 1e-3, "drainage area not close enough")
  expect(max(abs(catchment_area$totda - catchment_area$nhdptotda)) < 1e-2, "drainage area not close enough")
})

test_that("arbolate sum works", {
  source(system.file("extdata", "walker_data.R", package = "nhdplusTools"))
  catchment_length <- prepare_nhdplus(walker_flowline, 0, 0,
                                      purge_non_dendritic = FALSE, warn = FALSE)

  catchment_length <- select(walker_flowline, COMID) %>%
    left_join(catchment_length, by = "COMID") %>%
    select(ID = COMID, toID = toCOMID, length = LENGTHKM)

  arb_sum <- calculate_arbolate_sum(catchment_length)

  catchment_length$arb_sum <- arb_sum
  catchment_length$nhd_arb_sum <- walker_flowline$ArbolateSu

  expect(mean(abs(catchment_length$arb_sum - catchment_length$nhd_arb_sum)) < 1e-3, "arbolate sum not close enough")
  expect(max(abs(catchment_length$arb_sum - catchment_length$nhd_arb_sum)) < 1e-2, "arbolate sum not close enough")
})

test_that("get_pfaf", {
  suppressMessages(
  source(system.file("extdata/nhdplushr_data.R", package = "nhdplusTools")))
  hr_flowline <- align_nhdplus_names(hr_data$NHDFlowline)

  suppressWarnings(
  fl <- prepare_nhdplus(hr_flowline, 0, 0, purge_non_dendritic = FALSE, warn = FALSE))

  fl <- select(hr_flowline, COMID, AreaSqKM) %>%
    left_join(fl, by = "COMID") %>%
    st_sf() %>%
    select(ID = COMID, toID = toCOMID, area = AreaSqKM)

  fl$nameID = " "
  fl$totda <- calculate_total_drainage_area(sf::st_set_geometry(fl, NULL))
  fl <- left_join(fl, get_levelpaths(dplyr::rename(sf::st_set_geometry(fl, NULL),
                                                   weight = totda)), by = "ID")

  pfaf <- get_pfaf(fl, max_level = 2)

  expect_equal(pfaf[pfaf$ID == 15000500028335,	], dplyr::tibble(ID = 15000500028335,
                                                         pf_level_1 = 5, pf_level_2 = 51))

  pfaf <- get_pfaf(fl, max_level = 4)

  expect_true(all(!is.na(c(pfaf$pf_level_1, pfaf$pf_level_4))))

  fl <- left_join(fl, pfaf, by = "ID")

  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500061836], 611)

  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500028338], 591)
  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500050711], 592)
  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500028337], 593)

  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500072804], 151)
  expect_equal(pfaf$pf_level_4[pfaf$ID == 15000500072804], 1511)

  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500084318], 161)
  expect_equal(pfaf$pf_level_3[pfaf$ID == 15000500028332], 181)

  source(system.file("extdata", "walker_data.R", package = "nhdplusTools"))

  fl <- prepare_nhdplus(walker_flowline, 0, 0, purge_non_dendritic = FALSE, warn = FALSE)

  fl <- select(walker_flowline, COMID, AreaSqKM) %>%
    left_join(fl, by = "COMID") %>%
    st_sf() %>%
    select(ID = COMID, toID = toCOMID, area = AreaSqKM)

  fl$nameID = ""
  fl$totda <- calculate_total_drainage_area(sf::st_set_geometry(fl, NULL))
  fl <- left_join(fl, get_levelpaths(dplyr::rename(sf::st_set_geometry(fl, NULL),
                                                   weight = totda)), by = "ID")

  pfaf <- get_pfaf(fl, max_level = 2)

  expect_equal(nrow(pfaf), 57)
})

test_that("get_terminal", {
  suppressMessages(
    source(system.file("extdata/nhdplushr_data.R", package = "nhdplusTools")))
  hr_flowline <- align_nhdplus_names(hr_data$NHDFlowline)

  suppressWarnings(
    fl <- prepare_nhdplus(hr_flowline, 0, 0, purge_non_dendritic = FALSE, warn = FALSE) %>%
      select(ID = COMID, toID = toCOMID))

  outlet <- fl$ID[which(!fl$toID %in% fl$ID)]
  fl$toID[which(!fl$toID %in% fl$ID)] <- 0
  terminal <- get_terminal(fl, outlet)

  expect_equal(names(terminal), c("terminalID", "ID"))
  expect_true(is.numeric(terminal$ID))
  expect_equal(nrow(terminal), nrow(fl))

  source(system.file("extdata", "walker_data.R", package = "nhdplusTools"))

  fl <- prepare_nhdplus(walker_flowline, 0, 0, purge_non_dendritic = FALSE, warn = FALSE) %>%
    select(ID = COMID, toID = toCOMID)

  outlet <- fl$ID[which(!fl$toID %in% fl$ID)]
  fl$toID[which(!fl$toID %in% fl$ID)] <- 0

  terminal <- get_terminal(fl, outlet)

  expect_equal(nrow(terminal), nrow(fl))
  expect_true(is.integer(terminal$ID))
})

test_that("get_terminal", {
  source(system.file("extdata", "walker_data.R", package = "nhdplusTools"))

  fl <- prepare_nhdplus(walker_flowline, 0, 0, purge_non_dendritic = FALSE, warn = FALSE) %>%
    select(ID = COMID, toID = toCOMID, length = LENGTHKM)

  suppressWarnings(pl <- get_pathlength(fl))

  expect_equal(nrow(fl), nrow(pl))

  pl <- left_join(select(walker_flowline,
                         COMID, Pathlength),
                  pl,
                  by = c("COMID" = "ID"))

  expect_equal(pl$pathlength, pl$Pathlength)
})
