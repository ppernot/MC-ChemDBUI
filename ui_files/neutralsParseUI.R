tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Parse",.noWS = "outside"),
      tabsetPanel(
        tabPanel(
          title = 'Select',
          fluidRow(
            column(
              8,
              shiny::textInput(
                "neutralsReacSel",
                "Species filter"
              )
            ),
            column(
              4,
              actionButton(
                "neutralsReacSelInit",
                "Reset"
              ),
              tags$style(
                type='text/css',
                "#neutralsReacSelInit { width:100%; margin-top: 30px;}"
              )
            )
          ),
          radioButtons(
            "neutralsReacSelKind",
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
          uiOutput("selNeutralsReac"),
          hr(),
          actionButton(
            "neutralsParseSave",
            "Apply changes",
            icon = icon('save',verify_fa = FALSE)
          )
        ),
        tabPanel(
          title = 'Plot',
          fluidRow(
            column(
              8,
              selectInput(
                'neutralsSimulateSize',
                label    = '# MC samples',
                choices  = seq(100,1000,by = 100),
                selected = 500
              )
            ),
            column(
              4,
              actionButton(
                "neutralsSimulateBtn",
                label = "Go !",
                icon  = icon('gear',verify_fa = FALSE),
                class = "btn-primary"
              ),
              tags$style(
                type='text/css',
                "#neutralsSimulateBtn { width:100%; margin-top: 30px;}"
              )
            )
          ),
          sliderInput(
            "neutralsTempRangePlot",
            label = "Temp. range [K]",
            min   = 10, 
            max   = 600,
            value = c(50,350),
            step  = 50,
            round = TRUE
          ),
          numericInput(
            "M0Plot",
            label = "Ref. Log10(density [cm^-3])",
            value = 18,
            min = 5, max = 22, step = 1,
            width = '150px'
          ),
          hr(),
          sliderInput(
            "densRangePlot",
            label = "Log10(density [cm^-3]) range",
            min   = 5, 
            max   = 22,
            value = c(8,20),
            step  = 1,
            round = TRUE
          ),
          numericInput(
            "T0Plot",
            label = "Ref. Temp. [K]",
            value = 150,
            min = 10, max = 600, step = 10,
            width = '150px'
          )
        )      
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            'Rate params',
            fluidRow(
              column(
                8,
                uiOutput("neutralsRateMask")
              ),
              column(
                4,
                plotOutput("plotNeutralsRate",height = plotHeight)
              )
            )
          ),
          # tabPanel(
          #   'Biblio',
          #   uiOutput("neutralsBiblio")
          # ),
          tabPanel(
            'Help',
            fluidRow(
              column(
                6,
                HTML(
                  "<h4>Neutrals-Parse</h4> Visualize and edit data for individual
                  reactions.
                  <ul>
                    <li> <strong>Search</strong>: enter a species name to 
                  filter reactions. Reinitialize with Reset button.
                    <li> <strong>Reactions</strong>: list of reactions in DB, 
                  possibly filtered. The up and down arrows enable to go through
                  the list step by step.
                  </ul>
                  <h4>Simulation</h4> Generate random samples to
                  build graphs for rate constants and branching ratios.
                  <br>
                  <h4>Apply changes</h4> click to apply the 
                  changes made to the reaction's data. To save to disk, go to 
                  the Edit page. 
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
