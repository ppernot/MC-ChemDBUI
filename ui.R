function(request) {
  source_ui <- function(...) {
    source(
      file.path("ui_files", ...),
      local = TRUE
    )$value
  }

  navbarPage(
    "MC-ChemDB",
    theme = shinythemes::shinytheme(
      c("cosmo", "cerulean", "spacelab", "yeti")[3]
    ),
    navbarMenu(
      "Neutrals",
      tabPanel(
        title = "Edit",
        source_ui("neutralsEdit.R")
      ),
      tabPanel(
        title = "Parse",
        source_ui("neutralsParse.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("neutralsSample.R")
      )
    ),
    tabPanel(
      title = "Ions",
      source_ui("ions.R")
    ),
    tabPanel(
      title = "Photo-processes",
      source_ui("photo.R")
    ),
    tabPanel(
      title = "About",
      source_ui("about.R")
    )
  )
}
