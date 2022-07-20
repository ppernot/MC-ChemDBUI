tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Edit",.noWS = "outside"),
      fluidRow(
        column(
          6,
          uiOutput("selNeuOrigVersion"),
        ),
        column(
          6,
          textInput(
            "neuCopyVersion",
            "Target DB Version:",
            value = NULL,
            placeholder = "v_X.X"
          )
        )
      ),
      uiOutput("selNeuFile"),
      actionButton(
        "neuSave",
        "Save",
        icon = icon('save',verify_fa = FALSE)
      ),
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
  
