tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Edit/View",.noWS = "outside"),
      tabsetPanel(
        tabPanel(
          title = 'Select',
          br(),
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
          fluidRow(
            column(
              8,
              shiny::selectizeInput(
                "neutralsReaction",
                "Reactions",
                choices  = NULL,
                options = list(maxOptions = maxOptions)
              )
            ),
            column(
              4,
              fluidRow(
                column(
                  6,
                  actionButton(
                    "neutralsMinus",
                    "",
                    icon = icon('angle-down',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#neutralsMinus { width:100%; margin-top: 30px;}"
                  )
                ),
                column(
                  6,
                  actionButton(
                    "neutralsPlus",
                    "",
                    icon = icon('angle-up',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#neutralsPlus { width:100%; margin-top: 30px;}"
                  )
                )
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              6,
              checkboxInput(
                "neutralsParseComment",
                "Comment reaction",
                value = FALSE
              )
            ),
            column(
              6,
              actionButton(
                "neutralsParseSave",
                "Apply changes",
                icon = icon('save',verify_fa = FALSE)
              )
            )
          )
        ),
        tabPanel(
          title = 'Plot',
          br(),
          selectInput(
            'neutralsSimulateSize',
            label    = '# MC samples',
            choices  = seq(100,1000,by = 100),
            selected = 500,
            width = '200px'
          ),
          hr(),
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
          tabPanel(
            'Biblio',
            uiOutput("neutralsBiblio")
          ),
          tabPanel(
            'Help',
            fluidRow(
              column(
                6,
                HTML(
                  "<h4>Neutrals-Edit/View</h4> Visualize and edit data for individual
                  reactions.
                  <ul>
                    <li> <strong>Search</strong>: enter a species name to 
                  filter reactions. Reinitialize with Reset button.
                    <li> <strong>Reactions</strong>: list of reactions in DB, 
                  possibly filtered. The up and down arrows enable to go through
                  the list step by step.
                    <li> <strong>Comment reaction</strong>: when selected, 
                  the reaction will be commented out and not parsed.
                    <li> <strong>Apply changes</strong>: click to apply the 
                  changes made to the reaction's data. To save to disk, go to 
                  the Load page. 
                  </ul>
                  <h4>Plot</h4> Parameters to to build graphs for rate laws,
                  using random nMC samples.
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
