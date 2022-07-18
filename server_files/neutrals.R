neutralsVersion = shiny::reactiveVal()
neutralsFile    = shiny::reactiveVal()
neutralsReacs   = shiny::reactiveVal()

output$selNeuVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "neuVersion",
      "Neutrals DB Version:",
      rev(
        list.dirs(
          path=neutralsSource, 
          full.names = FALSE, 
          recursive = FALSE)
      )
    )
  )
})

output$selNeuFile = shiny::renderUI({
  req(neutralsVersion())
  list(
    shiny::selectInput(
      "neuFile",
      "Neutrals DB File:",
      c("Choose a file..." = "",
        list.files(
          path=file.path(neutralsSource,neutralsVersion()), 
          full.names = FALSE, 
          recursive = FALSE)
      )
    )
  )
})

shiny::observe({
  neutralsVersion(input$neuVersion)
})

shiny::observe({
  neutralsFile(input$neuFile)
})

shiny::observe({
  req(neutralsVersion())
  req(neutralsFile())
  neutralsReacs(
    readLines(
      con <- file(file.path(
        neutralsSource,
        neutralsVersion(),
        neutralsFile()
      ))
    )
  )
  close(con)
})

output$neuHeader = shiny::renderUI({
  list(
    h4(file.path(neutralsVersion(),neutralsFile()))
  )
})
  
shiny::observe({
  shinyAce::updateAceEditor(
    session, "ace", 
    value = paste(neutralsReacs(),collapse = '\n')
  )
})
