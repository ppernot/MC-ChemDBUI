tabPanel(
  title = "Load",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Load",.noWS = "outside"),
      fluidRow(
        column(
          6,
          uiOutput("selNeuOrigVersion"),
        ),
        column(
          6,
          textInput(
            "neutralsCopyVersion",
            "Target DB Version:",
            value = NULL,
            placeholder = "v_X.X"
          )
        )
      ),
      # uiOutput("selNeuFile"),
      fluidRow(
        column(
          6,
          actionButton(
            "neutralsEditSave",
            "Save",
            icon = icon('save',verify_fa = FALSE)
          )
        ),
        column(
          6,
          actionButton(
            "neutralsEditRestore",
            "Restore",
            icon = icon('trash-undo',verify_fa = FALSE)
          )
        )
      ),
      hr(),
      h4("Statistics"),
      verbatimTextOutput(
        "neuDBStats",
        placeholder = TRUE
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            "Database",
            uiOutput("neutralsHeaderDB"),
            aceEditor(
              outputId = "aceNeutralsDB",
              cursorId = "cursor",
              selectionId = "selection",
              height = "600px",
              fontSize = 16,
              showPrintMargin = FALSE,
              placeholder = "Please select a DB version..."
            )
          ),
          tabPanel(
            "Release Notes",
            uiOutput("neutralsHeaderRN"),
            aceEditor(
              outputId = "aceNeutralsRN",
              cursorId = "cursor",
              selectionId = "selection",
              height = "600px",
              fontSize = 16,
              showPrintMargin = FALSE,
              placeholder = "Please select a DB version..."
            )
          )
        )
      )
    )
  )
)
  
