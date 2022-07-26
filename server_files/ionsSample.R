output$selIonsVersionSample = shiny::renderUI({
  list(
    shiny::selectInput(
      "ionsVersionSample",
      "Source DB Version:",
      rev(
        list.dirs(
          path=ionsSource, 
          full.names = FALSE, 
          recursive  = FALSE)
      )
    )
  )
})

observeEvent(
  input$ionsSampleBtn,
  {
    fp = file.path(ionsSource,input$ionsVersionSample,'Data')
    listDirs = list.files(
      path       = fp,
      full.names = FALSE, 
      recursive  = FALSE
    )
    dataFiles = paste0(fp,'/',listDirs,'/data.csv')
    finf1     = file.info(dataFiles, extra_cols = FALSE)
    fpTmp     = file.path(ionsTmp,input$ionsVersionSample,'Reactions')
    tmpFiles  = paste0(fpTmp,'/',listDirs,'/Samples/run_0000.csv')
    finf2     = file.info(tmpFiles, extra_cols = FALSE)
    selByMod  = which(difftime(finf1[,"mtime"],finf2[,"mtime"]) > 0)
    listDirs  = listDirs[selByMod]
    print(listDirs)
  })
