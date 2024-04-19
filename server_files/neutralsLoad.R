neutralsEditOrigVersion = shiny::reactiveVal()
neutralsEditCopyVersion = shiny::reactiveVal()
neutralsEditDBFile      = shiny::reactiveVal('neutralsDB.csv')
neutralsEditDBFileBack  = shiny::reactiveVal('neutralsDB.bck')
neutralsEditDBText      = shiny::reactiveVal()
neutralsEditRNFile      = shiny::reactiveVal('ReleaseNotes.txt')
neutralsEditRNFileBack  = shiny::reactiveVal('ReleaseNotes.bck')
neutralsEditRNText      = shiny::reactiveVal()
neutralsDB              = shiny::reactiveVal()
neutralsDBStats         = shiny::reactiveVal()

# Select ####
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

# Load ####
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
  data[1] = paste0(data[1],',"ID" ')
  for (i in 2:length(data))
    data[i] = paste0(data[i],',',i-1)

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

# Parse ####
makeTag = function(x) 
  trimws(paste0(x[1],': ',x[2],' -> ',x[3]))

shiny::observe({
  req(input$aceNeutralsDB)
  neutralsDBStats(NULL) # Reset stats while processing
  
  # Get all data in (do not use comment.char here because it mixes up with id.)
  # data = try(
  #   read.table(
  #     header = TRUE, 
  #     text = input$aceNeutralsDB, 
  #     sep=';',
  #     quote = "\"",
  #     comment.char = ""),
  #   silent = TRUE
  # )
  data = try(
    data.table::fread(
      text = input$aceNeutralsDB 
    ),
    silent = TRUE
  )
  if(class(data) == "try-error") {
    id = shiny::showNotification(
      strong(paste0('Cannot parse DB: incorrect format! Msg:',data)),
      closeButton = TRUE,
      duration = NULL,
      type = 'error'
    )
    return(NULL)
  }
  
  # Filter out comment lines
  sel = substr(data$R1,1,1) != '#'
  data = data[sel,]
  
  data[['REACTANTS']] = apply(data[,paste0("R",1:3)], 1, 
                              function(x) paste0(x[!is.na(x) & x != ""],
                                                 collapse = ' + '))
  data[['PRODUCTS']]  = apply(data[,paste0("P",1:5)], 1, 
                              function(x) paste0(x[!is.na(x) & x != ""],
                                                 collapse = ' + '))
  # Build unique tag
  data[['TAG']]       = apply(data[,c("ID","REACTANTS","PRODUCTS")], 1, 
                              makeTag)
  
  # Check that reactions are mass balanced
  reacList = prodList = c()
  for(i in seq_along(data$REACTANTS)) {
    reactants = getSpecies(data$REACTANTS[i])
    reacList  = c(reacList, reactants)
    products  = getSpecies(data$PRODUCTS[i])
    prodList  = c(prodList, products)
    msg = checkBalance(reactants,products)
    if(!is.null(msg))
      id = shiny::showNotification(
        strong(paste0(msg,', at line ',data$ID[i])),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
    req(is.null(msg))
  }
  
  stats = list(
    nbReac = nrow(data),
    nbReactants = length(unique(reacList)),
    nbProducts = length(unique(prodList)),
    nbSpecies = length(unique(c(reacList,prodList)))
  )
  neutralsDBStats(stats)
  
  neutralsDB(data)
})
output$neuDBStats = shiny::renderText({
  req(neutralsDBStats())
  paste0(
    'Nb. reactions = ',neutralsDBStats()$nbReac,'\n',
    'Nb. species   = ',neutralsDBStats()$nbSpecies,'\n',
    'incl.\n',
    '    reactants = ',neutralsDBStats()$nbReactants,'\n',
    '    products  = ',neutralsDBStats()$nbProducts
  )
})

# Save ####
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
        # Remove last column (ID)
        S = unlist(strsplit(x, ','))
        paste0(S[-length(S)], collapse = ',')
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

# Restore ####
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
