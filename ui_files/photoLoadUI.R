tabPanel(
  title = "Files",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("PhotoProcs - Files",.noWS = "outside"),
      fluidRow(
        column(
          6,
          uiOutput("selPhotoEditOrigVersion"),
        ),
        column(
          6,
          textInput(
            "photoEditCopyVersion",
            "Target DB Version:",
            value = NULL,
            placeholder = "v_X.X"
          )
        )
      ),
      fluidRow(
        column(
          6,
          actionButton(
            "photoEditSave",
            "Save",
            icon = icon('save',verify_fa = FALSE)
          )
        ),
        column(
          6,
          actionButton(
            "photoEditRestore",
            "Restore",
            icon = icon('trash-undo',verify_fa = FALSE)
          )
        )
      ),
      hr(),
      h4("Statistics"),
      verbatimTextOutput(
        "photoDBStats",
        placeholder = TRUE
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            "Scheme",
            uiOutput("photoHeaderDB"),
            aceEditor(
              outputId = "acePhotoDB",
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
            uiOutput("photoHeaderRN"),
            aceEditor(
              outputId = "acePhotoRN",
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
