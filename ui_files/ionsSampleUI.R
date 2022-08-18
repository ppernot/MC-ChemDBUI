tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Ions - Sample",.noWS = "outside"),
      uiOutput("selIonsVersionSample"),
      fluidRow(
        column(
          6,
          checkboxInput(
            'ionsSampleUpdate',
            label = 'Update',
            value = TRUE
          )
        ),
        column(
          6,
          checkboxInput(
            'ionsSampleCheck',
            label = 'Check only',
            value = FALSE
          )
        )
      ),
      fluidRow(
        column(
          7,
          selectInput(
            'ionsSampleSize',
            label    = '# MC samples',
            choices  = c(0,10,100,500,1000),
            selected = 500
          )
        ),
        column(
          5,
          actionButton(
            "ionsSampleBtn",
            label = "Sample !",
            icon  = icon('gear',verify_fa = FALSE),
            class = "btn-primary"
          ),
          tags$style(
            type='text/css',
            "#ionsSampleBtn { width:100%; margin-top: 30px;}"
          )#,
          # actionButton(
          #   "ionsCnvrtBtn",
          #   label = "Convert !",
          #   icon  = icon('gear',verify_fa = FALSE),
          #   class = "btn-primary"
          # ),
          # tags$style(
          #   type='text/css',
          #   "#ionsCnvrtBtn { width:100%; margin-top: 30px;}"
          # )
        )
      )
    ),
    mainPanel(
      width = mainWidth,
      wellPanel(
        tabsetPanel(
          tabPanel(
            title = "Statistics",
            verbatimTextOutput("ionsStats",)
          ),
          tabPanel(
            title = "Help",
            HTML(
              "<h4>Sample</h4> 
          Generate MC samples. <BR>
          For efficiency, this is done in two steps:
          <OL>
          <LI> generate intermediate samples for all reactions, stored in temporary files.
          <LI> gather and collate intermediate samples to ChemDBPublic.
          </OL>
          This enables DB updates without regenerating all samples, which might be  
          processor intensive.
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
