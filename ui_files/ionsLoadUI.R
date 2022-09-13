tabPanel(
  title = "Files",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Files",.noWS = "outside"),
      fluidRow(
        column(
          6,
          uiOutput("selIonsEditOrigVersion"),
        ),
        column(
          6,
          textInput(
            "ionsEditCopyVersion",
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
            "ionsEditSave",
            "Save",
            icon = icon('save',verify_fa = FALSE)
          )
        ),
        column(
          6,
          actionButton(
            "ionsEditRestore",
            "Restore",
            icon = icon('trash-undo',verify_fa = FALSE)
          )
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            "Database",
            uiOutput("ionsHeaderDB"),
            aceEditor(
              outputId = "aceIonsDB",
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
            uiOutput("ionsHeaderRN"),
            aceEditor(
              outputId = "aceIonsRN",
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
