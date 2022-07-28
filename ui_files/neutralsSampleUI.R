tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Sample"),
      tabsetPanel(
        tabPanel(
          "Generate",
          br(),
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
                label = "Go !",
                icon  = icon('gear',verify_fa = FALSE),
                class = "btn-primary"
              ),
              tags$style(
                type='text/css',
                "#neutralsSampleBtn { width:100%; margin-top: 20px;}"
              )
            )
          ),
          wellPanel(
            h4("About",.noWS = "before"),
            HTML("
                 Generate random samples from DB prescriptions. 
                 This may take some time...
                 ")
          )
        ),
        tabPanel(
          "Plot",
          br(),
          textInput(
            "reacNbPlot",
            label = "Reaction id",
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
          ),
          wellPanel(
            h4("About",.noWS = "before"),
            HTML("
                 Plot T-dependent and density-dependent realizations
                 of the generated samples. Controls define the range 
                 of plots and the reference density and temperature
                 used for the cuts in the (T,density) plane.<br>
                 <i>Rq</i>: need to <b>Parse</b> first... 
                 ")
          )  
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      plotOutput("plotRate")
    )
  )
)
