photoEditOrigVersion = shiny::reactiveVal()
photoEditCopyVersion = shiny::reactiveVal()
photoEditDBFile      = shiny::reactiveVal('photoDB.csv')
photoEditDBFileBack  = shiny::reactiveVal('photoDB.bck')
photoEditDBText      = shiny::reactiveVal()
photoEditRNFile      = shiny::reactiveVal('ReleaseNotes.txt')
photoEditRNFileBack  = shiny::reactiveVal('ReleaseNotes.bck')
photoEditRNText      = shiny::reactiveVal()
photoDB              = shiny::reactiveVal()
photoDBStats         = shiny::reactiveVal()

output$selPhotoEditOrigVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "photoEditOrigVersion",
      "Source DB Version:",
      rev(
        list.dirs(
          path = photoSource, 
          full.names = FALSE, 
          recursive  = FALSE)
      )
    )
  )
})

shiny::observe({
  photoEditOrigVersion(input$photoEditOrigVersion)
})

shiny::observe({
  photoEditCopyVersion(input$photoEditCopyVersion)
})

# Load ####
shiny::observe({
  req(photoEditOrigVersion())
  req(photoEditDBFile())
    
  data = readLines(
    con <- file(file.path(
      photoSource,
      photoEditOrigVersion(),
      photoEditDBFile()
    ))
  )
  close(con)
  
  # Add a unique id to each line for editing
  data[1] = paste0('ID;',data[1])
  for (i in 2:length(data))
    data[i] = paste0(i-1,';',data[i])
  
  photoEditDBText(data)
  
})

shiny::observe({
  req(photoEditOrigVersion())
  req(photoEditRNFile())
  photoEditRNText(
    readLines(
      con <- file(file.path(
        photoSource,
        photoEditOrigVersion(),
        photoEditRNFile()
      ))
    )
  )
  close(con)
})

output$photoHeaderDB = shiny::renderUI({
  list(
    strong(file.path(photoEditOrigVersion(),photoEditDBFile()))
  )
})

output$photoHeaderRN = shiny::renderUI({
  list(
    strong(file.path(photoEditOrigVersion(),photoEditRNFile()))
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "acePhotoDB", 
    value = paste(photoEditDBText(),collapse = '\n')
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "acePhotoRN", 
    value = paste(photoEditRNText(),collapse = '\n')
  )
})

# Parse ####
shiny::observe({
  req(input$acePhotoDB)
  
  # Get all data in (do not use comment.char here because it mixes up with id.)
  data = try(
    read.table(
      header = TRUE, 
      text = input$acePhotoDB, 
      sep = ';',
      quote = "\"",
      comment.char = ""),
    silent = TRUE
  )
  
  if(class(data)=="try-error") {
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
  
  data[['REACTANTS']] = apply(data[,paste0("R",1:2)], 1, 
                              function(x) paste0(x[!is.na(x) & x != ""],
                                                 collapse = ' + '))
  data[['PRODUCTS']]  = apply(data[,paste0("P",1:4)], 1, 
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
  photoDBStats(stats)
  
  photoDB(data)
})

output$photoDBStats = shiny::renderText({
  req(photoDBStats())
  paste0(
    'Nb. reactions = ',photoDBStats()$nbReac,'\n',
    'Nb. species   = ',photoDBStats()$nbSpecies,'\n',
    'incl.\n',
    '    reactants = ',photoDBStats()$nbReactants,'\n',
    '    products  = ',photoDBStats()$nbProducts
  )
})

# Save ####
shiny::observeEvent(
  input$photoEditSave,
  {
    req(photoEditCopyVersion())

    photoOrigDir = file.path(photoSource,photoEditOrigVersion())
    photoCopyDir = file.path(photoSource,photoEditCopyVersion())
    if(!dir.exists(photoCopyDir)) {
      dir.create(photoCopyDir)
      id = shiny::showNotification(
        h4(paste0('Created new version: ', photoEditCopyVersion())),
        closeButton = FALSE,
        duration = 5
      )
    }

    if(photoOrigDir == photoCopyDir) {
      # Create backup files
      file.copy(
        from = file.path(photoOrigDir,photoEditDBFile()),
        to   = file.path(photoOrigDir,photoEditDBFileBack())
      )
      file.copy(
        from = file.path(photoOrigDir,photoEditRNFile()),
        to   = file.path(photoOrigDir,photoEditRNFileBack())
      )
    }
    
    # Save DB and RN files to target version
    # - for DB, need to remove ID before saving
    dataEditor = photoEditDBText()
    data = sapply(
      dataEditor, 
      function(x) {
        i = regexpr(";",x)
        substr(x,i+1,nchar(x))
      }
    )
    writeLines(
      data,
      con = file.path(photoCopyDir,photoEditDBFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', photoEditDBFile())),
      closeButton = FALSE,
      duration = 5
    )

    data = isolate(input$acePhotoRN)
    writeLines(
      data,
      con = file.path(photoCopyDir,photoEditRNFile())
    )
    id = shiny::showNotification(
      h4(paste0('Saved file: ', photoEditRNFile())),
      closeButton = FALSE,
      duration = 5
    )
    
  }
)

# Restore ####
shiny::observeEvent(
  input$photoEditRestore,
  {
    req(photoEditCopyVersion())
    req(photoEditOrigVersion() == photoEditCopyVersion())
    photoOrigDir = file.path(photoSource,photoEditOrigVersion())
    req(file.exists(file.path(photoOrigDir,photoEditDBFileBack())))
    req(file.exists(file.path(photoOrigDir,photoEditRNFileBack())))
    
    shiny::showModal(shiny::modalDialog(
      title = "This will erase all changes since last Save.",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("photoEditRestoreOK", "OK")
      )
    ))
  }
)

shiny::observeEvent(
  input$photoEditRestoreOK,
  {
    req(photoEditCopyVersion())
    req(photoEditOrigVersion() == photoEditCopyVersion())
    
    photoOrigDir = file.path(photoSource,photoEditOrigVersion())
    
    photoEditDBText(
      readLines(
        con <- file(file.path(photoOrigDir,photoEditDBFileBack()))
      )
    )
    close(con)
    photoEditRNText(
      readLines(
        con <- file(file.path(photoOrigDir,photoEditRNFileBack()))
      )
    )
    close(con)
    removeModal()
  }
)
