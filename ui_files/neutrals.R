sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    uiOutput("selNeuVersion"),
    uiOutput("selNeuFile")
  ),
  mainPanel(
    width = mainWidth,
    uiOutput("neuHeader"),
    aceEditor(
      outputId = "ace",
      # to access content of `selectionId` in server.R use `ace_selection`
      # i.e., the outputId is prepended to the selectionId for use
      # with Shiny modules
      selectionId = "selection",
      height = "600px",
      fontSize = 14,
      showPrintMargin = FALSE,
      placeholder = "Please select a DB version and a DB file..."
    )
  )
)
