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
                "photoBRSampleSort",
                label = "Sort samples",
                value = TRUE
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
            #   column(
            #     6,
            #     checkboxInput(
            #       'photoSampleUpdate',
            #       label = 'Update',
            #       value = TRUE
            #     )
            #   ),
            #   column(
            #     6,
            #     checkboxInput(
            #       'photoSampleCheck',
            #       label = 'Check only',
            #       value = FALSE
            #     )
            #   )
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
          tabPanel(
            title = "Statistics",
            verbatimTextOutput("photoStats",)
          ),
          tabPanel(
            title = "Help",
            HTML(
              "<h4>Sample</h4> 
          Generate MC samples. <BR>
          <h4>Options</h4>
          <UL>
          <LI> <B>Update</B>: if checked, only reactions recently modified or with missing
          samples are sampled. If unchecked, all reactions are sampled, which might 
          take some time...
          <LI> <B>Check only</B>: if checked, the data are tested for consistency, 
          but samples are not generated.
          <LI> <B># MC samples</B>: number of samples, default (and recommended minimum
          for production) is 500. Smaller values are used for testing. Zero (0) 
          corresponds to a single draw with nominal values of parameters. It is always 
          provided as the first sample in file run_0000.csv.
          </UL>"
            )
          )
        )  
      )
    )
  )
)
