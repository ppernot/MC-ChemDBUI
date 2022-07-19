tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Sample"),
      selectInput(
        'sampleSize',
        label    = '# MC samples (0: nominal)',
        choices  = c(seq(0,100,by=10),seq(200,1000,by = 100)),
        selected = 0
      ),
      actionButton(
        "neutralsSampleBtn",
        label = "Generate",
        icon  = icon('gear',verify_fa = FALSE)
      ),
      wellPanel(
        textInput(
          "reacNbPlot",
          label = NULL,
          value = NULL,
          placeholder = "Enter reac. number"
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      plotOutput("plotRate")
    )
  )
)
