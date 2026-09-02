if(identical(Sys.getenv("NOT_CRAN"), "true")) {
  library("testthat")

  test_check("nhdplusTools")
}
