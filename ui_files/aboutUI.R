tabPanel(
  title = "About",
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tagList(
        h5(paste0("Version : ",version)),
        # h5(a(href="https://github.com/ppernot/MC-ChemDBUI","How to cite")),
        h5(a(href="https://github.com/ppernot/MC-ChemDBUI","code@github")),
        h5(a(href="https://github.com/ppernot/MC-ChemDBUI/issues",
             "Bug report, Feature request")),
        hr( style="border-color: #666;"),
        h5("Author : P. Pernot"),
        h5("Affiliation : CNRS"),
        h5(a(href="https://github.com/ppernot","Contact")),
        # h5(a(href="https://ppernot.github.io/UncVal","User's Manual")),
      )
    ),
    mainPanel()
  )
)
