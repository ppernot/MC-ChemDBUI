neutralsOrigVersion = shiny::reactiveVal()
neutralsCopyVersion = shiny::reactiveVal()
neutralsFile        = shiny::reactiveVal()
neutralsReacs       = shiny::reactiveVal()

output$selNeuOrigVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "neuOrigVersion",
      "Source DB Version:",
      rev(
        list.dirs(
          path=neutralsSource, 
          full.names = FALSE, 
          recursive  = FALSE)
      )
    )
  )
})

output$selNeuFile = shiny::renderUI({
  req(neutralsOrigVersion())
  list(
    shiny::selectInput(
      "neuFile",
      "Neutrals DB File:",
      c("Choose a file..." = "",
        list.files(
          path=file.path(neutralsSource,neutralsOrigVersion()), 
          full.names = FALSE, 
          recursive = FALSE)
      )
    )
  )
})

shiny::observe({
  neutralsOrigVersion(input$neuOrigVersion)
})

shiny::observe({
  neutralsCopyVersion(input$neuCopyVersion)
})

shiny::observe({
  neutralsFile(input$neuFile)
})

shiny::observe({
  req(neutralsOrigVersion())
  req(neutralsFile())
  neutralsReacs(
    readLines(
      con <- file(file.path(
        neutralsSource,
        neutralsOrigVersion(),
        neutralsFile()
      ))
    )
  )
  close(con)
})

output$neuHeader = shiny::renderUI({
  list(
    h4(file.path(neutralsOrigVersion(),neutralsFile()))
  )
})
  
shiny::observe({
  shinyAce::updateAceEditor(
    session, "ace", 
    value = paste(neutralsReacs(),collapse = '\n')
  )
})

output$checkSpecies <- renderText({ 
  req(input$ace_cursor)
  sp = input$ace_selection
  compo = get.atoms(sp)
  names(compo) = elements
  mass  = massFormula(compo)
  paste0(
    "Selection: \"", sp, "\"\n",
    "Mass= ", mass
  )
})

shiny::observeEvent(
  input$neuSave,
  {
    req(neutralsCopyVersion())
    req(neutralsFile())
   
    # Make copy of source directory
    neutralsOrigDir = file.path(neutralsSource,neutralsOrigVersion())
    neutralsCopyDir = file.path(neutralsSource,neutralsCopyVersion())
    if(!dir.exists(neutralsCopyDir)) {
      dir.create(neutralsCopyDir)
      files = list.files(path = neutralsOrigDir)
      for(file in files)
        file.copy(
          from = file.path(neutralsOrigDir,file),
          to   = file.path(neutralsCopyDir,file)
        )
      id = shiny::showNotification(
        h4(paste0('Created new version: ', neutralsCopyVersion())),
        closeButton = FALSE,
        duration = 5
      )
    }
    
    # Save modified file to target version
    data = isolate(input$ace)
    writeLines(
      data,
      con = file.path(neutralsCopyDir,neutralsFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', neutralsFile())),
      closeButton = FALSE,
      duration = 5
    )
    
  }
)
