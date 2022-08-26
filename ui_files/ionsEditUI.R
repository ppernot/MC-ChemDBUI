tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Edit/View",.noWS = "outside"),
      tabsetPanel(
        tabPanel(
          'Select',
          br(),
          fluidRow(
            column(
              8,
              shiny::textInput(
                "ionsReacSel",
                "Species filter"
              )
            ),
            column(
              4,
              actionButton(
                "ionsReacSelInit",
                "Reset"
              ),
              tags$style(
                type='text/css',
                "#ionsReacSelInit { width:100%; margin-top: 30px;}"
              )
            )
          ),
          radioButtons(
            "ionsReacSelKind",
            label = "",
            choices = c(
              "Reactant" = "Reactant",
              "Product" = "Product",
              "Both" = "Both"
            ),
            selected = "Both",
            inline = TRUE
          ),
          br(),
          uiOutput("selIonsReac"),
          hr(),
          actionButton(
            "ionsParseSave",
            "Apply changes",
            icon = icon('save',verify_fa = FALSE)
          )
        ),
        tabPanel(
          'Plot',
          selectInput(
            'ionsSimulateSize',
            label    = '# MC samples',
            choices  = seq(100,1000,by = 100),
            selected = 500,
            width = '200px'
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
        )
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
                4,
                uiOutput("ionsBRMask")
              ),
              column(
                8,
                tabsetPanel(
                  tabPanel(
                    'Parallel',
                    plotOutput("plotIonsBRSample",height = plotHeight)
                  ),
                  tabPanel(
                    'Tree',
                    plotOutput("plotIonsBRTree",height = plotHeight)
                  )
                )
              )
            )
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
                  "<h4>Ions-Edit/View</h4> Visualize and edit data for individual
                  reactions.
                  <ul>
                    <li> <strong>Search</strong>: enter a species name to 
                  filter reactions. Reinitialize with Reset button.
                    <li> <strong>Reactions</strong>: list of reactions in DB, 
                  possibly filtered. The up and down arrows enable to go through
                  the list step by step.
                    <li> <strong>Apply changes</strong>: click to apply the 
                  changes made to the reaction's data. To save to disk, go to 
                  the Load page.
                  </ul>
                  <h4>Plot</h4> Generate random samples and
                  build graphs for rate constants and branching ratios.
                  <br>
                   
                  "
                )
              )
            )
          )
        )
      )
    )
  )
)
