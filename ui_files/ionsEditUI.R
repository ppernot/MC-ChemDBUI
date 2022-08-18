tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Edit",.noWS = "outside"),
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
      uiOutput("selIonsEditFile"),
      actionButton(
        "ionsEditSave",
        "Save",
        icon = icon('save',verify_fa = FALSE)
      )  
    ),
    mainPanel(
      width = mainWidth,
      uiOutput("ionsHeader"),
      aceEditor(
        outputId = "aceIons",
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
