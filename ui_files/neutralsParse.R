tabPanel(
  title = "Parse",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Parse"),
      actionButton(
        "neutralsParseBtn",
        label = "Parse",
        icon  = icon('gear',verify_fa = FALSE)
      )
    ),
    mainPanel(
      width = mainWidth,
      tabsetPanel(
        tabPanel(
          "Reactions",
          br(),
          fluidRow(
            column(
              3,
              textInput(
                "targetSpecies",
                label = NULL,
                value = NA,
                placeholder = "Filter by species"
              )
            ),
            column(
              6,
              radioButtons(
                "targetSpeciesKind",
                label = "",
                choices = c(
                  "Reactant" = "Reactant",
                  "Product" = "Product",
                  "Both" = "Both"
                ),
                selected = "Both",
                inline = TRUE
              )
            )
          ),
          DT::DTOutput("tabScheme")
        ),
        tabPanel(
          title = "Checks",
          # h4("Species with no mass"),
          # verbatimTextOutput("spNoMass"),
          h4("Species per mass"),
          DT::DTOutput("massScheme")
        )
      )
    )
  )
)
