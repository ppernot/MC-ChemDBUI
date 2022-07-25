tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Edit",.noWS = "outside"),
      fluidRow(
        column(
          6,
          uiOutput("selIonsOrigVersion"),
        ),
        column(
          6,
          textInput(
            "ionsCopyVersion",
            "Target DB Version:",
            value = NULL,
            placeholder = "v_X.X"
          )
        )
      ),
      uiOutput("selIonsFile"),
      actionButton(
        "ionsSave",
        "Save",
        icon = icon('save',verify_fa = FALSE)
      ),
      hr(),
      h4('Simulation'),
      fluidRow(
        column(
          8,
          selectInput(
            'ionsSampleSize',
            label    = '# MC samples (0: nominal)',
            choices  = c(seq(0,100,by=10),seq(200,1000,by = 100)),
            selected = 100
          )
        ),
        column(
          4,
          actionButton(
            "ionsSimulateBtn",
            label = "Go !",
            icon  = icon('gear',verify_fa = FALSE),
            class = "btn-primary"
          ),
          tags$style(
            type='text/css',
            "#ionsSimulateBtn { width:100%; margin-top: 20px;}"
          )
        )
      ),
      sliderInput(
        "ionsTempRangePlot",
        label = "Temp. range [K]",
        min   = 100, 
        max   = 1000,
        value = c(100,1000),
        step  = 100,
        round = TRUE
      ),
      br(),
      wellPanel(
        h4("About",.noWS = "before"),
        HTML("<B>Edit</B> Choose a DB version and a DB file.
             <br>
             <B>Simulation</B> Random samples are generated to
             build graphs for rate constants and branching 
             ratios (Graphs panel)
             ")
      )  
    ),
    mainPanel(
      width = mainWidth,
      uiOutput("ionsHeader"),
      tabsetPanel(
        tabPanel(
          'Mask',
          uiOutput("ionsMask")
        ),
        tabPanel(
          'Graphs',
          tabsetPanel(
            tabPanel(
              'Rate',
              plotOutput("plotIonsParsSample")
            ),
            tabPanel(
              'BR-Sample',
              textAreaInput(
                "ionsStringDist",
                width = '600px',
                cols  = 120,
                label = 'Distribution'
              ),
              plotOutput("plotIonsBRSample")
            ),
            tabPanel(
              'BR-Tree',
              plotOutput("plotIonsBRTree")
            )
          )
        ),
        tabPanel(
          'Text',
          aceEditor(
            outputId = "aceIons",
            cursorId = "cursor",
            selectionId = "selection",
            height = "600px",
            fontSize = 16,
            showPrintMargin = FALSE,
            useSoftTabs     = FALSE, # Do NOT replace tabs by spaces
            showInvisibles  = TRUE,
            placeholder = "Please select a DB version and a DB file..."
          )
        )
      )
    )
  )
)
