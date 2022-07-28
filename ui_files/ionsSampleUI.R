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
            label    = '# MC samples (0:nominal)',
            choices  = c(0,10,seq(100,1000,by = 100)),
            selected = 0
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
            "#ionsSampleBtn { width:100%; margin-top: 20px;}"
          )
        )
      ),
      br(),
      wellPanel(
        h4("About",.noWS = "before"),
        HTML(
          "<B>Sample</B> Generate MC samples. For efficiency, this is done in
          two steps:
          <OL>
          <LI> Generate intermediate samples for all reactions, stored in temporary files.
          <LI> Gather and collate intermediate samples to ChemDBPublic
          </OL>
          This enables to update the DB without regenerating all samples.
          <UL>
          <LI> <B>Update</B>: if checked, only reactions recently modified or with missing
          samples are sampled. If unchecked, all reactions are sampled, which might 
          take some time...
          <LI> <B>Check only</B>: if checked, the data are tested for consistency, but samples are
          not generated.
          </UL>
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
