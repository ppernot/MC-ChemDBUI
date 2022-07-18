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
      ),
      hr(),
      textInput(
        "targetSpecies",
        label = NULL,
        value = NA,
        placeholder = "Select species"
      )
    ),
    mainPanel(
      width = mainWidth,
      tabsetPanel(
        tabPanel(
          title = "Reactions", 
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
