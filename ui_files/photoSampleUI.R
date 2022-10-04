tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("PhotoProcs - Sample",.noWS = "outside"),
      br(),
      # fluidRow(
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
      # ),
      fluidRow(
        column(
          7,
          selectInput(
            'photoSampleSize',
            label    = '# MC samples',
            choices  = c(0,10,100,500,1000),
            selected = 500
          )
        ),
        column(
          5,
          actionButton(
            "photoSampleBtn",
            label = "Sample !",
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
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
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
