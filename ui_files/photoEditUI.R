tabPanel(
  title = "Edit",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("PhotoProcs - Edit/View",.noWS = "outside"),
      tabsetPanel(
        tabPanel(
          'Select',
          br(),
          fluidRow(
            column(
              8,
              shiny::textInput(
                "photoReacSel",
                "Species filter"
              )
            ),
            column(
              4,
              actionButton(
                "photoReacSelInit",
                "Reset"
              ),
              tags$style(
                type='text/css',
                "#photoReacSelInit { width:100%; margin-top: 30px;}"
              )
            )
          ),
          br(),
          fluidRow(
            column(
              8,
              shiny::selectizeInput(
                "photoReaction",
                "Reactions",
                choices = NULL,
                options = list(maxOptions = maxOptions)
              )
            ),
            column(
              4,
              fluidRow(
                column(
                  6,
                  actionButton(
                    "photoMinus", "",
                    icon = icon('angle-down',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#photoMinus { width:100%; margin-top: 30px;}"
                  )
                ),
                column(
                  6,
                  actionButton(
                    "photoPlus", "",
                    icon = icon('angle-up',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#photoPlus { width:100%; margin-top: 30px;}"
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
                "photoParseComment",
                "Comment reaction",
                value = FALSE
              )
            ),
            column(
              6,
              actionButton(
                "photoParseSave",
                "Apply changes",
                icon = icon('save',verify_fa = FALSE)
              )
            )
          )
        ),
        tabPanel(
          'Plot',
          selectInput(
            'photoSimulateSize',
            label    = '# MC samples',
            choices  = seq(100,1000,by = 100),
            selected = 100,
            width = '200px'
          ),
          selectInput(
            "photoXSReso",
            "Resolution (nm):",
            photoXSResolutions,
            selected = 1
          ),
          sliderInput(
            "photoWLPlotRange",
            label = "Wavelength [nm]",
            min   = 50, 
            max   = 350,
            value = c(50,250),
            step  =  10,
            round = TRUE
          ),
          checkboxInput(
            "photoBRSort",
            label = "Sort samples",
            value = TRUE
          )
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            'Cross-section',
            fluidRow(
              column(
                5,
                uiOutput("photoXSMaskUI")
              ),
              column(
                7,
                plotOutput("plotPhotoXSSample",height = plotHeight)
              )
            )
          ),
          tabPanel(
            'BRs',
            fluidRow(
              column(
                5,
                uiOutput("photoBRMaskUI")
              ),
              column(
                7,
                plotOutput("plotPhotoBRSample",height = plotHeight)
              )
            )
          ),
          tabPanel(
            'Biblio',
            uiOutput("photoBiblio")
          ),
          tabPanel(
            'Help',
            fluidRow(
              column(
                6,
                HTML(
                  "<h4>Photo-Edit/View</h4> Visualize and edit data for individual
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
