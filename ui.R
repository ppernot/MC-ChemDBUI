function(request) {
  navbarPage(
    "MC-ChemDB",
    theme = bslib::bs_theme(
      version = 5, 
      bootswatch = c("united","sketchy")[1]
    ),
    navbarMenu(
      "Neutrals",
      tabPanel(
        title = "Edit",
        source_ui("neutralsEditUI.R")
      ),
      tabPanel(
        title = "Parse",
        source_ui("neutralsParseUI.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("neutralsSampleUI.R")
      ),
      tabPanel(
        title = "Report",
        source_ui("neutralsReportUI.R")
      )
    ),
    navbarMenu(
      "Ions",
      tabPanel(
        title = "Edit",
        source_ui("ionsEditUI.R")
      ),
      tabPanel(
        title = "Parse",
        source_ui("ionsParseUI.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("ionsSampleUI.R")
      ),
      tabPanel(
        title = "Report",
        source_ui("ionsReportUI.R")
      )
    ),
    tabPanel(
      title = "Photo-processes",
      source_ui("photoUI.R")
    ),
    tabPanel(
      title = "About",
      source_ui("aboutUI.R")
    )
  )
}
