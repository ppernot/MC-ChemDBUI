neutralsEditOrigVersion = shiny::reactiveVal()
neutralsEditCopyVersion = shiny::reactiveVal()
neutralsEditDBFile      = shiny::reactiveVal('neutralsDB.csv')
neutralsEditDBFileBack  = shiny::reactiveVal('neutralsDB.bck')
neutralsEditDBText      = shiny::reactiveVal()
neutralsEditRNFile      = shiny::reactiveVal('ReleaseNotes.txt')
neutralsEditRNFileBack  = shiny::reactiveVal('ReleaseNotes.bck')
neutralsEditRNText      = shiny::reactiveVal()

output$selNeuOrigVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "neutralsOrigVersion",
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

shiny::observe({
  neutralsEditOrigVersion(input$neutralsOrigVersion)
})

shiny::observe({
  neutralsEditCopyVersion(input$neutralsCopyVersion)
})

shiny::observe({
  req(neutralsEditOrigVersion())
  req(neutralsEditDBFile())
  
  data = readLines(
    con <- file(file.path(
      neutralsSource,
      neutralsEditOrigVersion(),
      neutralsEditDBFile()
    ))
  )
  close(con)
  
  # Add a unique id to each line for editing
  data[1] = paste0('ID;',data[1])
  for (i in 2:length(data))
    data[i] = paste0(i-1,';',data[i])

  neutralsEditDBText(data)
})

shiny::observe({
  req(neutralsEditOrigVersion())
  req(neutralsEditRNFile())
  neutralsEditRNText(
    readLines(
      con <- file(file.path(
        neutralsSource,
        neutralsEditOrigVersion(),
        neutralsEditRNFile()
      ))
    )
  )
  close(con)
})

output$neutralsHeaderDB = shiny::renderUI({
  list(
    strong(file.path(neutralsEditOrigVersion(),neutralsEditDBFile()))
  )
})

output$neutralsHeaderRN = shiny::renderUI({
  list(
    strong(file.path(neutralsEditOrigVersion(),neutralsEditRNFile()))
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceNeutralsDB", 
    value = paste(neutralsEditDBText(),collapse = '\n')
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceNeutralsRN", 
    value = paste(neutralsEditRNText(),collapse = '\n')
  )
})

shiny::observeEvent(
  input$neutralsEditSave,
  {
    req(neutralsEditCopyVersion())
    
    neutralsOrigDir = file.path(neutralsSource,neutralsEditOrigVersion())
    neutralsCopyDir = file.path(neutralsSource,neutralsEditCopyVersion())
    if(!dir.exists(neutralsCopyDir)) {
      dir.create(neutralsCopyDir)
      id = shiny::showNotification(
        h4(paste0('Created new version: ', neutralsEditCopyVersion())),
        closeButton = FALSE,
        duration = 5
      )
    }
    
    if(neutralsOrigDir == neutralsCopyDir) {
      # Create backup files
      file.copy(
        from = file.path(neutralsOrigDir,neutralsEditDBFile()),
        to   = file.path(neutralsOrigDir,neutralsEditDBFileBack())
      )
      file.copy(
        from = file.path(neutralsOrigDir,neutralsEditRNFile()),
        to   = file.path(neutralsOrigDir,neutralsEditRNFileBack())
      )
    }
    
    # Save DB and RN files to target version
    # For DB, need to remove ID before saving
    dataEditor = neutralsEditDBText()
    data = sapply(
      dataEditor, 
      function(x) {
        i = regexpr(";",x)
        substr(x,i+1,nchar(x))
      }
    )
    writeLines(
      data,
      con = file.path(neutralsCopyDir,neutralsEditDBFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', neutralsEditDBFile())),
      closeButton = FALSE,
      duration = 5
    )
    
    data = isolate(input$aceNeutralsRN)
    writeLines(
      data,
      con = file.path(neutralsCopyDir,neutralsEditRNFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', neutralsEditRNFile())),
      closeButton = FALSE,
      duration = 5
    )
    
  }
)

shiny::observeEvent(
  input$neutralsEditRestore,
  {
    req(neutralsEditCopyVersion())
    req(neutralsEditOrigVersion() == neutralsEditCopyVersion())
    neutralsOrigDir = file.path(neutralsSource,neutralsEditOrigVersion())
    req(file.exists(file.path(neutralsOrigDir,neutralsEditDBFileBack())))
    req(file.exists(file.path(neutralsOrigDir,neutralsEditRNFileBack())))
    
    shiny::showModal(shiny::modalDialog(
      title = "This will erase all changes since last Save.",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("neutralsEditRestoreOK", "OK")
      )
    ))
  }
)

shiny::observeEvent(
  input$neutralsEditRestoreOK,
  {
    req(neutralsEditCopyVersion())
    req(neutralsEditOrigVersion() == neutralsEditCopyVersion())
    
    neutralsOrigDir = file.path(neutralsSource,neutralsEditOrigVersion())
    
    neutralsEditDBText(
      readLines(
        con <- file(file.path(neutralsOrigDir,neutralsEditDBFileBack()))
      )
    )
    close(con)
    neutralsEditRNText(
      readLines(
        con <- file(file.path(neutralsOrigDir,neutralsEditRNFileBack()))
      )
    )
    close(con)
    removeModal()
  }
)
