ionsEditOrigVersion = shiny::reactiveVal()
ionsEditCopyVersion = shiny::reactiveVal()
ionsEditDBFile      = shiny::reactiveVal('ionsDB.csv')
ionsEditDBFileBack  = shiny::reactiveVal('ionsDB.bck')
ionsEditDBText      = shiny::reactiveVal()
ionsEditRNFile      = shiny::reactiveVal('ReleaseNotes.txt')
ionsEditRNFileBack  = shiny::reactiveVal('ReleaseNotes.bck')
ionsEditRNText      = shiny::reactiveVal()

output$selIonsEditOrigVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "ionsEditOrigVersion",
      "Source DB Version:",
      rev(
        list.dirs(
          path = ionsSource, 
          full.names = FALSE, 
          recursive  = FALSE)
      )
    )
  )
})

shiny::observe({
  ionsEditOrigVersion(input$ionsEditOrigVersion)
})

shiny::observe({
  ionsEditCopyVersion(input$ionsEditCopyVersion)
})

shiny::observe({
  req(ionsEditOrigVersion())
  req(ionsEditDBFile())
  ionsEditDBText(
    readLines(
      con <- file(file.path(
        ionsSource,
        ionsEditOrigVersion(),
        ionsEditDBFile()
      ))
    )
  )
  close(con)
})

shiny::observe({
  req(ionsEditOrigVersion())
  req(ionsEditRNFile())
  ionsEditRNText(
    readLines(
      con <- file(file.path(
        ionsSource,
        ionsEditOrigVersion(),
        ionsEditRNFile()
      ))
    )
  )
  close(con)
})

output$ionsHeaderDB = shiny::renderUI({
  list(
    strong(file.path(ionsEditOrigVersion(),ionsEditDBFile()))
  )
})

output$ionsHeaderRN = shiny::renderUI({
  list(
    strong(file.path(ionsEditOrigVersion(),ionsEditRNFile()))
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceIonsDB", 
    value = paste(ionsEditDBText(),collapse = '\n')
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceIonsRN", 
    value = paste(ionsEditRNText(),collapse = '\n')
  )
})

shiny::observeEvent(
  input$ionsEditSave,
  {
    req(ionsEditCopyVersion())

    ionsOrigDir = file.path(ionsSource,ionsEditOrigVersion())
    ionsCopyDir = file.path(ionsSource,ionsEditCopyVersion())
    if(!dir.exists(ionsCopyDir)) {
      dir.create(ionsCopyDir)
      id = shiny::showNotification(
        h4(paste0('Created new version: ', ionsEditCopyVersion())),
        closeButton = FALSE,
        duration = 5
      )
    }

    if(ionsOrigDir == ionsCopyDir) {
      # Create backup files
      file.copy(
        from = file.path(ionsOrigDir,ionsEditDBFile()),
        to   = file.path(ionsOrigDir,ionsEditDBFileBack())
      )
      file.copy(
        from = file.path(ionsOrigDir,ionsEditRNFile()),
        to   = file.path(ionsOrigDir,ionsEditRNFileBack())
      )
    }
    
    # Save DB and RN files to target version
    data = isolate(input$aceIonsDB)
    writeLines(
      data,
      con = file.path(ionsCopyDir,ionsEditDBFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', ionsEditDBFile())),
      closeButton = FALSE,
      duration = 5
    )

    data = isolate(input$aceIonsRN)
    writeLines(
      data,
      con = file.path(ionsCopyDir,ionsEditRNFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', ionsEditRNFile())),
      closeButton = FALSE,
      duration = 5
    )
    
  }
)

shiny::observeEvent(
  input$ionsEditRestore,
  {
    req(ionsEditCopyVersion())
    req(ionsEditOrigVersion() == ionsEditCopyVersion())
    ionsOrigDir = file.path(ionsSource,ionsEditOrigVersion())
    req(file.exists(file.path(ionsOrigDir,ionsEditDBFileBack())))
    req(file.exists(file.path(ionsOrigDir,ionsEditRNFileBack())))
    
    shiny::showModal(shiny::modalDialog(
      title = "This will erase all changes since last Save.",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ionsEditRestoreOK", "OK")
      )
    ))
  }
)

shiny::observeEvent(
  input$ionsEditRestoreOK,
  {
    req(ionsEditCopyVersion())
    req(ionsEditOrigVersion() == ionsEditCopyVersion())
    
    ionsOrigDir = file.path(ionsSource,ionsEditOrigVersion())
    
    ionsEditDBText(
      readLines(
        con <- file(file.path(ionsOrigDir,ionsEditDBFileBack()))
      )
    )
    close(con)
    ionsEditRNText(
      readLines(
        con <- file(file.path(ionsOrigDir,ionsEditRNFileBack()))
      )
    )
    close(con)
    removeModal()
  }
)
