tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Sample"),
      fluidRow(
        column(
          8,
          selectInput(
            'sampleSize',
            label    = '# MC samples (0: nominal)',
            choices  = c(seq(0,100,by=10),seq(200,1000,by = 100)),
            selected = 0
          )
        ),
        column(
          4,
          actionButton(
            "neutralsSampleBtn",
            label = "Generate",
            icon  = icon('gear',verify_fa = FALSE),
            class = "btn-primary"
          ),
          tags$style(
            type='text/css',
            "#neutralsSampleBtn { width:100%; margin-top: 20px;}"
          )
        )
      ),
      br(),
      wellPanel(
        h4("Plots"),
        textInput(
          "reacNbPlot",
          label = NULL,
          value = NULL,
          placeholder = "Enter reac. number",
          width = '150px'
        ),
        sliderInput(
          "tempRangePlot",
          label = "Temp. range [K]",
          min   = 10, 
          max   = 600,
          value = c(50,350),
          step  = 50,
          round = TRUE
        ),
        textInput(
          "M0Plot",
          label = "Ref. Log10(density [cm^-3])",
          value = 18,
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
        textInput(
          "T0Plot",
          label = "Ref. Temp. [K]",
          value = 150,
          width = '150px'
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      plotOutput("plotRate")
    )
  )
)
