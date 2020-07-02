#' @title shiny app LC-N2G
#' @export
li <- function() {
  appDir <- system.file("LC-N2G", package = "LCN2G")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

