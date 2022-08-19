ionsEditOrigVersion = shiny::reactiveVal()
ionsEditCopyVersion = shiny::reactiveVal()
ionsEditDBFile      = shiny::reactiveVal('ionsDB.csv')
ionsEditDBText      = shiny::reactiveVal()
ionsEditRNFile      = shiny::reactiveVal('ReleaseNotes.txt')
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

# output$selIonsEditFile = shiny::renderUI({
#   req(ionsEditOrigVersion())
#   list(
#     shiny::selectInput(
#       "ionsEditFile",
#       "File to edit:",
#       c("Choose a file..." = "",
#         list.files(
#           path=file.path(ionsSource,ionsEditOrigVersion()), 
#           full.names = FALSE, 
#           recursive = FALSE)
#       )
#     )
#   )
# })

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
    req(ionsEditFile())

    # Make copy of source directory
    ionsOrigDir = file.path(ionsSource,ionsEditOrigVersion())
    ionsCopyDir = file.path(ionsSource,ionsEditCopyVersion())
    if(!dir.exists(ionsCopyDir)) {
      dir.create(ionsCopyDir)
      # files = list.files(path = ionsOrigDir)
      # for(file in files)
      #   file.copy(
      #     from = file.path(ionsOrigDir,file),
      #     to   = file.path(ionsCopyDir,file)
      #   )
      id = shiny::showNotification(
        h4(paste0('Created new version: ', ionsEditCopyVersion())),
        closeButton = FALSE,
        duration = 5
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
