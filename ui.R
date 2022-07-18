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
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    tabPanel(
      title = "Neutrals",
      source_ui("neutrals.R")
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
