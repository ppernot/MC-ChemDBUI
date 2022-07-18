tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Edit",.noWS = "outside"),
      uiOutput("selNeuVersion"),
      uiOutput("selNeuFile"),
      h5("Check Species Mass"),
      verbatimTextOutput("checkSpecies", placeholder = TRUE)
    ),
    mainPanel(
      width = mainWidth,
      uiOutput("neuHeader"),
      aceEditor(
        outputId = "ace",
        cursorId = "cursor",
        selectionId = "selection",
        height = "600px",
        fontSize = 16,
        showPrintMargin = FALSE,
        placeholder = "Please select a DB version and a DB file..."
      )
    )
  )
)
  
