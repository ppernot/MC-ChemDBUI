tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Parse",.noWS = "outside"),
      actionButton(
        "ionsParseSave",
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
                plotOutput("plotIonsBRSample",height = plotHeight)
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
            'Help',
            fluidRow(
              column(
                6,
                HTML(
                  "<h4>Parse</h4> Choose a DB version and a reaction.
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
