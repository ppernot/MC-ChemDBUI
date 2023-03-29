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
          ),
          hr(),
          checkboxInput(
            "photoEditAdvanced",
            label = "Advanced options",
            value = FALSE
          ),
          conditionalPanel(
            condition = "input.photoEditAdvanced",
            conditionalPanel(
              condition = "!input.photoEditGP_Fit",
              fluidRow(
                column(
                  4,
                  checkboxInput(
                    "useDirg",
                    label = "Dirg (vs. Diri)",
                    value = TRUE
                  )
                ),
                column(
                  4,
                  conditionalPanel(
                    condition = "!input.useDirg",
                    checkboxInput(
                      "newDiri",
                      label = "New Diri",
                      value = TRUE
                    )
                  ),
                  conditionalPanel(
                    condition = "input.useDirg",
                    checkboxInput(
                      "newDirg",
                      label = "New Dirg",
                      value = TRUE
                    )
                  )
                  
                  
                ),
                column(
                  4,
                  checkboxInput(
                    "flatTree",
                    label = "Flat tree",
                    value = TRUE
                  )
                ),
              ),
              fluidRow(
                column(
                  6,
                  checkboxInput(
                    "photoBRArrange",
                    label = "Arrange samples",
                    value = TRUE
                  )
                ),
                column(
                  6,
                  conditionalPanel(
                    condition = "input.photoBRArrange",
                    checkboxInput(
                      "photoBRUseRanks",
                      label = "Use ranks",
                      value = TRUE
                    )
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "!input.photoBRArrange",
              fluidRow(
                column(
                  6,
                  checkboxInput(
                    "photoEditGP_Fit",
                    label = "Gaussian Process"
                  ) 
                )
              ),
              conditionalPanel(
                condition = "input.photoEditGP_Fit",
                sliderInput(
                  "photoGPCorLen",
                  label = "GP corr. len. [nm]",
                  min   =   5, 
                  max   = 100,
                  value =  20,
                  step  =   5,
                  round = TRUE
                ),
                sliderInput(
                  "photoGPCorVar",
                  label = "GP Var.",
                  min   =   1, 
                  max   =  20,
                  value =   2,
                  step  =   1,
                  round = TRUE
                )
              )
            )
          )
        ),
        tabPanel(
          'Plot',
          br(),
          fluidRow(
            column(
              6,
              selectInput(
                'photoSimulateSize',
                label    = '# MC samples',
                choices  = c(10,seq(100,1000,by = 100)),
                selected = 10,
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
            min   = 10, 
            max   = 550,
            value = c(50,250),
            step  =  10,
            round = TRUE
          ),
          fluidRow(
            column(
              6,
              selectInput(
                'photoEditBRDisplay',
                label    = 'Display of Brs',
                choices  = c(
                  "All channels"   = 0,
                  "Neus vs Ions"   = 1,
                  "Sum-to-one"     = 2,
                  "Rel. uncert."    = 3,
                  "Mean rel. unc."   = 4
                )
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
                shinycssloaders::withSpinner(
                  plotOutput("plotPhotoXSSample",height = plotHeight)
                )
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
                shinycssloaders::withSpinner(
                  plotOutput("plotPhotoBRSample",height = plotHeight)
                )
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
                  <li> <strong>Advanced options</strong>
                      <ul> 
                        <li><strong>Dirg</strong> Use Dirg instead of Diri
                        <li><strong>New Dirg/Diri</strong> Use improved versions
                           of the Diri or Dirg distributions
                        <li><strong>Flat tree</strong> do not use a nested model
                           to separate ions from neutral channels (default: TRUE)
                        <li><strong>Arrange samples</strong> reorder the samples
                           to increase wavelength-wise correlation (default: TRUE)
                        <li><strong>Use ranks</strong> base reordering of samples
                           on their ranks instead of their values
                      </ul>
                  </ul>
                  <h5>Plot</h5>
                  <ul>
                    <li> <strong># MC samples</strong>: number of samples 
                         to plot
                    <li> <strong>Resolution</strong>: wavelength resolution
                    <li> <strong>Wavelength</strong>: defines the wavelength range 
                         of the plot 
                    
                    <li> <strong>log XS</strong>: use a log axis for 
                         cross-sections.
                   <li> <strong>Display of BRs</strong>:
                         <ul>
                            <li> All Channels: all channels
                            <li> Neus vs Ions: sum of neutral
                                 channels and  sum of ionic channels
                            <li> Sum-to-one: sum of all channels
                                 (should be 1 if sampling is OK)
                            <li> Rel. uncert: relative uncertainty 
                                 of all channels
                            <li> Mean rel. unc.: mean relative uncertainty 
                                 over ions, neutrals and all channels
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
