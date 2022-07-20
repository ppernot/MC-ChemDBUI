tabPanel(
  title = "Parse",
  sidebarLayout(
    sidebarPanel(
      width = sideWidth,
      h4("Neutrals - Parse"),
      fluidRow(
        column(
          3,offset = 8,
          actionButton(
            "neutralsParseBtn",
            label = "Parse !",
            icon  = icon('gear',verify_fa = FALSE),
            class = "btn-primary"
          )
        )
      ),
      br(),br(),
      wellPanel(
        h4("About",.noWS = "before"),
        HTML("Parse the data files to extract and check information
             (reactants, products, rate law type, parameters...).
             The reactions are checked for mass balance.<br>
             The reaction list might be filtered by species,
             either as reactant, product or both.<br>
             A complementary check is provided as a list of species
             sorted by mass. All isomers or states of a species should 
             fall at the same mass.
             ")
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
              2,
              textInput(
                "targetSpecies",
                label = NULL,
                value = NA,
                placeholder = "Filter by species"
              )
            ),
            column(
              8,
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
          h4("Species vs. mass"),
          DT::DTOutput("massScheme")
        )
      )
    )
  )
)
