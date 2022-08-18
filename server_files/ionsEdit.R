ionsEditOrigVersion = shiny::reactiveVal()
ionsEditCopyVersion = shiny::reactiveVal()
ionsEditFile        = shiny::reactiveVal()
ionsEditText        = shiny::reactiveVal()

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

output$selIonsEditFile = shiny::renderUI({
  req(ionsEditOrigVersion())
  list(
    shiny::selectInput(
      "ionsEditFile",
      "File to edit:",
      c("Choose a file..." = "",
        list.files(
          path=file.path(ionsSource,ionsEditOrigVersion()), 
          full.names = FALSE, 
          recursive = FALSE)
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
  ionsEditFile(input$ionsEditFile)
})

shiny::observe({
  req(ionsEditOrigVersion())
  req(ionsEditFile())
  ionsEditText(
    readLines(
      con <- file(file.path(
        ionsSource,
        ionsEditOrigVersion(),
        ionsEditFile()
      ))
    )
  )
  close(con)
})

output$ionsHeader = shiny::renderUI({
  list(
    h4(file.path(ionsEditOrigVersion(),ionsEditFile()))
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceIons", 
    value = paste(ionsEditText(),collapse = '\n')
  )
})


shiny::observeEvent(
  input$ionsEditSave,
  {
    req(ionsEditCopyVersion())
    req(ionsEditFile())

    # Make copy of source directory
    ionsOrigDir = file.path(ionsSource,ionsEditOrigVersion())
    ionsCopyDir = file.path(ionsSource,ionsEditCopyVersion())
    if(!dir.exists(ionsCopyDir)) {
      dir.create(ionsCopyDir)
      files = list.files(path = ionsOrigDir)
      for(file in files)
        file.copy(
          from = file.path(ionsOrigDir,file),
          to   = file.path(ionsCopyDir,file)
        )
      id = shiny::showNotification(
        h4(paste0('Created new version: ', ionsEditCopyVersion())),
        closeButton = FALSE,
        duration = 5
      )
    }

    # Save modified file to target version
    data = isolate(input$aceIons)
    writeLines(
      data,
      con = file.path(ionsCopyDir,ionsEditFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', ionsEditFile())),
      closeButton = FALSE,
      duration = 5
    )

  }
)
