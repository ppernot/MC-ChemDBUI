tabPanel(
  title = "Sample",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Sample"),
      selectInput(
        'sampleSize',
        label    = '# MC samples (0: nominal)',
        choices  = seq(0,1000,by=100),
        selected = 500
      ),
      actionButton(
        "neutralsSampleBtn",
        label = "Generate",
        icon  = icon('gear',verify_fa = FALSE)
      )
    ),
    mainPanel(
      width = mainWidth,
      # tabsetPanel(
      #   tabPanel(
      #     "Reactions",
      #     br(),
      #     fluidRow(
      #       column(
      #         3,
      #         textInput(
      #           "targetSpecies",
      #           label = NULL,
      #           value = NA,
      #           placeholder = "Filter by species"
      #         )
      #       ),
      #       column(
      #         6,
      #         radioButtons(
      #           "targetSpeciesKind",
      #           label = "",
      #           choices = c(
      #             "Reactant" = "Reactant",
      #             "Product" = "Product",
      #             "Both" = "Both"
      #           ),
      #           selected = "Both",
      #           inline = TRUE
      #         )
      #       )
      #     ),
      #     DT::DTOutput("tabScheme")
      #   ),
      #   tabPanel(
      #     title = "Checks",
      #     # h4("Species with no mass"),
      #     # verbatimTextOutput("spNoMass"),
      #     h4("Species per mass"),
      #     DT::DTOutput("massScheme")
      #   )
      # )
    )
  )
)
