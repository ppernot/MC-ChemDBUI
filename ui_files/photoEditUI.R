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
          fluidRow(
            column(
              6,
              selectInput(
                'photoSimulateSize',
                label    = '# MC samples',
                choices  = seq(100,1000,by = 100),
                selected = 100,
                width = '200px'
              )
            ),
            column(
              6,
              selectInput(
                "photoXSReso",
                "Resolution (nm):",
                photoXSResolutions,
                selected = 1
              )
            )
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
          fluidRow(
            column(
              6,
              checkboxInput(
                "photoBRSort",
                label = "Sort samples",
                value = TRUE
              )
            ),
            column(
              6,
              checkboxInput(
                "photoXSLog",
                label = "log XS",
                value = TRUE
              )
            )
          ),
          fluidRow(
            column(
              6,
              selectInput(
                'photoEditBRDisplay',
                label    = 'Display of Brs',
                choices  = c(
                  "All channels" = 0,
                  "Neus vs Ions" = 1,
                  "Sum-to-one"   = 2
                )
              )
            )
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
                12,
                HTML(
                  "<h4>PhotoProcs-Edit/View</h4> Visualize and edit data for individual
                  reactions.
                  <h5>Select</h5> 
                  <ul>
                    <li> <strong>Species filter</strong>: enter a species name to 
                  filter reactions. Reinitialize with Reset button.
                    <li> <strong>Reactions</strong>: list of reactions in DB, 
                  possibly filtered. The up and down arrows enable to go through
                  the list step by step.
                    <li> <strong>Comment reaction</strong>: inactivate this reaction
                  in the database.
                    <li> <strong>Apply changes</strong>: click to apply the 
                  changes made to the reaction's data. To save to disk, go to 
                  the Load page.
                  </ul>
                  <h5>Plot</h5>
                  <ul>
                    <li> <strong># MC samples</strong>: number of samples 
                         to plot
                    <li> <strong>Résolution</strong>: wavelength resolution
                    <li> <strong>Wavelength</strong>: defines the wavelength range 
                         of the plot 
                    <li> <strong>Sort samples</strong>: used to preserve as much 
                         as possible a wavelength-wise continuity of the random 
                         samples. Introduces a (unknown) level of systematic
                         uncertainty wrt the pure random uncertainty generated 
                         by sampling.
                    <li> <strong>log XS</strong>: use a log axis for 
                         cross-sections.
                    <li> <strong>Display of BRs</strong>:
                         <ul>
                            <li> All Channels: display all channels
                            <li> Neus vs Ions: display the sum of neutral
                                 channels and the sum of ionic channels
                            <li> Sum-to-one: display the sum of all channels
                                 (should be 1 if sampling is OK)
                         </ul>
                  </ul>
                  
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
