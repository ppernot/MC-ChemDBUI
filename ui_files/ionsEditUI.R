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
            'ionsSimulateSize',
            label    = '# MC samples',
            choices  = seq(100,1000,by = 100),
            selected = 500
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
            "#ionsSimulateBtn { width:100%; margin-top: 30px;}"
          )
        )
      ),
      sliderInput(
        "ionsTempRangePlot",
        label = "Temp. range [K]",
        min   = 10, 
        max   = 600,
        value = c(100,500),
        step  =  50,
        round = TRUE
      )  
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        # uiOutput("ionsHeader"),
        tabsetPanel(
          tabPanel(
            'Rate',
            fluidRow(
              column(
                4,
                uiOutput("ionsRateMask")
              ),
              column(
                8,
                plotOutput("plotIonsParsSample",height = plotHeight)
              )
            )
          ),
          tabPanel(
            'BRs',
            fluidRow(
              column(
                6,
                uiOutput("ionsBRMask")
              ),
              column(
                6,
                tabPanel(
                  'BR-Sample',
                  textAreaInput(
                    "ionsStringDist",
                    width = '400px',
                    cols  = 120,
                    label = 'Distribution'
                  ),
                  plotOutput("plotIonsBRSample",height = plotHeight)
                )
              )
              
            )
          ),
          tabPanel(
            'BR-Tree',
            plotOutput("plotIonsBRTree",height = plotHeight)
          ),
          tabPanel(
            'Biblio',
            uiOutput("ionsBiblio")
          ),
          tabPanel(
            'CSV',
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
          ),
          tabPanel(
            'Help',
            fluidRow(
              column(
                6,
                HTML(
                  "<h4>Edit</h4> Choose a DB version and a DB file.
                <br>
                <h4>Simulation</h4> Random samples are generated to
                build graphs for rate constants and branching ratios."
                )
              )
            )
          )
        )
      )
    )
  )
)
