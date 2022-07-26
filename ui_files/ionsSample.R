tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Sample",.noWS = "outside"),
      uiOutput("selIonsVersionSample"),
      fluidRow(
        column(
          7,
          selectInput(
            'ionsSampleSize',
            label    = '# MC samples (0:nominal)',
            choices  = c(0,10,seq(100,1000,by = 100)),
            selected = 0
          )
        ),
        column(
          5,
          actionButton(
            "ionsSampleBtn",
            label = "Go !",
            icon  = icon('gear',verify_fa = FALSE),
            class = "btn-primary"
          ),
          tags$style(
            type='text/css',
            "#ionsSampleBtn { width:100%; margin-top: 20px;}"
          )
        )
      ),
      br(),
      wellPanel(
        h4("About",.noWS = "before"),
        HTML("<B>Edit</B> Choose a DB version and a DB file.
             <br>
             <B>Simulation</B> Random samples are generated to
             build graphs for rate constants and branching 
             ratios.
             ")
      )  
    ),
    mainPanel(
      width = mainWidth,
      tabsetPanel(
      )
    )
  )
)
