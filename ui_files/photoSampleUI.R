tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("PhotoProcs - Sample",.noWS = "outside"),
      tabsetPanel(
        tabPanel(
          "Generate",
          br(),
          fluidRow(
            column(
              6,
              checkboxInput(
                'photoSampleCheck',
                label = 'Check only',
                value = FALSE
              )
            ),
            column(
              6,
              selectInput(
                "photoSampleReso",
                "Resolution (nm):",
                c("All",photoXSResolutions),
                selected = "1"
              )
            )
          ),
          fluidRow(
            column(
              7,
              selectInput(
                'photoSampleSize',
                label    = '# MC samples',
                choices  = c(0,10,100,500,1000),
                selected = 10# 500
              )
            ),
            column(
              5,
              actionButton(
                "photoSampleBtn",
                label = "Go !",
                icon  = icon('gear',verify_fa = FALSE),
                class = "btn-primary"
              ),
              tags$style(
                type='text/css',
                "#photoSampleBtn { width:100%; margin-top: 30px;}"
              )
            )
          ),
          hr(),
          checkboxInput(
            "photoSampleAdvanced",
            label = "Advanced options",
            value = FALSE
          ),
          conditionalPanel(
            condition = "input.photoSampleAdvanced",
            fluidRow(
              column(
                6,
                checkboxInput(
                  "GP_Fit",
                  label = "Gaussian Process",
                  value = FALSE
                ) 
              ),
              column(
                6,
                checkboxInput(
                  "photoBRSampleSort",
                  label = "Sort samples",
                  value = FALSE
                )
              )
            )
          )
        ),
        tabPanel(
          "Plot",
          br(),
          fluidRow(
            column(
              8,
              shiny::selectizeInput(
                "photoSampleReaction",
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
                    "photoSampleMinus", "",
                    icon = icon('angle-down',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#photoSampleMinus { width:100%; margin-top: 30px;}"
                  )
                ),
                column(
                  6,
                  actionButton(
                    "photoSamplePlus", "",
                    icon = icon('angle-up',verify_fa = FALSE)
                  ),
                  tags$style(
                    type='text/css',
                    "#photoSamplePlus { width:100%; margin-top: 30px;}"
                  )
                )
              )
            )
          ),
          fluidRow(
            column(
              6,
              selectInput(
                'photoSamplePlotSize',
                label    = '# MC samples',
                choices  = c(0,10,100,500,1000),
                selected = 10# 500
              )
            ),
            column(
              6,
              selectInput(
                "photoSampleXSReso",
                "Resolution (nm):",
                photoXSResolutions,
                selected = 1
              )
            )
          ),
          sliderInput(
            "photoSampleWLPlotRange",
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
              selectInput(
                'photoSampleBRDisplay',
                label    = 'Display of Brs',
                choices  = c(
                  "All channels" = 0,
                  "Neus vs Ions" = 1,
                  "Sum-to-one"   = 2
                )
              )
            ),
            column(
              6,
              checkboxInput(
                "photoSampleXSLog",
                label = "log XS",
                value = TRUE
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
            title = "Plots",
            fluidRow(
              column(
                6,
                plotOutput("plotSampleXS", height = plotHeight)
              ),
              column(
                6,
                plotOutput("plotSampleBR", height = plotHeight)
              )
            )
          ),
          # tabPanel(
          #   title = "Statistics",
          #   verbatimTextOutput("photoStats",)
          # ),
          tabPanel(
            'Help',
            fluidRow(
              column(
                12,
                HTML(
                  "<h4>PhotoProcs-Sample</h4> Generate MC samples in ChemDBPublic.
                  
                  <h5>Generate</h5> 
                  <ul>
                    <li> <strong>Check only</strong>: goes through the process
                         without saving to disk. Can be used to verify that
                         generation will go smoothly...
                    <li> <strong>Sort samples</strong>: used to preserve as much 
                         as possible a wavelength-wise continuity of the random 
                         samples. Introduces a (unknown) level of systematic
                         uncertainty wrt the pure random uncertainty generated 
                         by sampling.
                    <li> <strong>Résolution</strong>: wavelength resolution
                    <li> <strong># MC samples</strong>: number of samples 
                         to generate
                    <li> <strong>Go !</strong>: click to start generation
                         process. This may take some time...
                    
                  </ul>
                  <h5>Plot</h5>
                  Read XS and BR samples on disk and plot them. 
                  <ul>
                    <li> <strong>Reactions</strong>: list of availables reactions
                    <li> <strong># MC samples</strong>: number of samples 
                         to plot
                    <li> <strong>Résolution</strong>: wavelength resolution
                    <li> <strong>Wavelength</strong>: defines the wavelength range 
                         of the plot 
                    
                    <li> <strong>Display of BRs</strong>:
                         <ul>
                            <li> All Channels: display all channels
                            <li> Neus vs Ions: display the sum of neutral
                                 channels and the sum of ionic channels
                            <li> Sum-to-one: display the sum of all channels
                                 (should be 1 if sampling is OK)
                         </ul>
                    <li> <strong>log XS</strong>: use a log axis for 
                         cross-sections.
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
