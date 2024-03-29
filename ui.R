function(request) {
  navbarPage(
    strong(paste0("MC-ChemDB ",version)),
    theme = bslib::bs_theme(
      version = 5, 
      bootswatch = c("united","sketchy")[1]
    ),
    
    navbarMenu(
      title = "PhotoProcs",
      tabPanel(
        title = "Files",
        source_ui("photoLoadUI.R")
      ),
      tabPanel(
        title = "Edit",
        source_ui("photoEditUI.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("photoSampleUI.R")
      )
    ),
    
    navbarMenu(
      "Neutrals",
      tabPanel(
        title = "Files",
        source_ui("neutralsLoadUI.R")
      ),
      tabPanel(
        title = "Edit",
        source_ui("neutralsEditUI.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("neutralsSampleUI.R")
      )
    ),
    
    navbarMenu(
      "Ions",
      tabPanel(
        title = "Files",
        source_ui("ionsLoadUI.R")
      ),
      tabPanel(
        title = "Edit",
        source_ui("ionsEditUI.R")
      ),
      tabPanel(
        title = "Sample",
        source_ui("ionsSampleUI.R")
      )
    ),
   
    tabPanel(
      title = "About",
      source_ui("aboutUI.R")
    )
  )
}
