# Sample ####
observeEvent(
  input$neutralsSampleBtn,
  {
    req(neutralsDB())
    
    # Output directory
    version = neutralsEditOrigVersion()
    if(!is.null(neutralsEditCopyVersion()) &
       neutralsEditCopyVersion() != "")
      version = neutralsEditCopyVersion()
    neutralsSampleDir = paste0(neutralsPublic,'_',version)
    if(!dir.exists(neutralsSampleDir))
      dir.create(neutralsSampleDir)
      
    nMC = as.numeric(input$sampleSize)
    
    shiny::withProgress(
      message = 'Sampling...', 
      {
        for (i in 0:nMC) {
          incProgress(1/nMC, detail = if(i==0) 'Nominal' else paste(i,'/',nMC))
          
          dbOut = data.frame(
            R1 = 'x', R2 = 'x', R3 = 'x',
            P1 = 'x', P2 = 'x', P3 = 'x', P4 = 'x',
            A1 = NA, B1 = NA, C1 = NA, F1 = NA, G1 = NA,
            A2 = NA, B2 = NA, C2 = NA, F2 = NA, G2 = NA,
            A3 = NA, B3 = NA, C3 = NA, F3 = NA, G3 = NA,
            FC = NA,
            TYPE = 'x'
          )
          
          for (m in 1:nrow(neutralsDB())) {

            pars = rep(0,length(neutralsRateParKwdList))
            names(pars) = neutralsRateParKwdList
            for (kwd in neutralsRateParKwdList)
              pars[kwd] = as.numeric(neutralsDB()[[kwd]][m])
            
            type = neutralsDB()$TYPE[m]
            
            sampleRateParams = oneSampleNeutralsPars(i, pars, type)
            names(sampleRateParams) = neutralsRateParKwdList
            
            dbOut[m,] = c(
              R1 = neutralsDB()$R1[m], 
              R2 = neutralsDB()$R2[m], 
              R3 = neutralsDB()$R3[m],
              P1 = neutralsDB()$P1[m], 
              P2 = neutralsDB()$P2[m], 
              P3 = neutralsDB()$P3[m], 
              P4 = neutralsDB()$P4[m],
              sampleRateParams,
              TYPE = type
            )
            
          }
          write.table(
            dbOut,
            file=file.path(
              neutralsSampleDir,
              paste0('run_',sprintf('%04i',i),'.csv')
            ),
            quote=TRUE, sep=';', na=" ",
            row.names=FALSE,col.names=FALSE)
        }
      })
  })

# Plot ####
observe({
  req(neutralsDB())
  tag = neutralsDB()$TAG
  nums = 1:length(tag)  
  names(nums) = tag
  shiny::updateSelectizeInput(
    session  = session,
    inputId  = "reacNbPlot", 
    choices  = nums,
    selected = "0",
    server   = TRUE
  )
})

output$plotRate = renderPlot({
  req(input$reacNbPlot)
  iReac = as.numeric(input$reacNbPlot)

  req(neutralsDB())
  tag     = neutralsDB()$TAG[[iReac]]
  params0 = rep(0,length(neutralsRateParKwdList))
  names(params0) = neutralsRateParKwdList
  for (kwd in neutralsRateParKwdList)
    params0[kwd] = as.numeric(neutralsDB()[[kwd]][iReac])
  typ0    = neutralsDB()$TYPE[[iReac]]
  note    = neutralsDB()$COMMENTS[[iReac]]

  legText = paste0(tag, '\n', 'Rate law: ', typ0)
  legText = switch(
    typ0,
    kooij    = paste0(
      legText,
      '\n',
      'Parameters: ',
      paste0(params0[1:5],
             collapse = ' / '),
      '\n'
    ),
    assocMD  = paste0(
      legText,
      '\n',
      'Parameters Fc : ', params0[16],'\n',
      'Parameters k0 : ',
      paste0(params0[1:5], collapse = ' / '),
      '\n',
      'Parameters kInf : ',
      paste0(params0[6:10], collapse = ' / '),
      '\n',
      'Parameters kr : ',
      paste0(params0[11:15], collapse = ' / '),
      '\n'
    ),
    assocVV  = paste0(
      legText,
      '\n',
      'Parameters Fc : ', params0[16],'\n',
      'Parameters kInf : ',
      paste0(params0[1:5], collapse = ' / '),
      '\n',
      'Parameters k0 : ',
      paste0(params0[6:10], collapse = ' / '),
      '\n',
      'Parameters kR : ',
      paste0(params0[11:15], collapse = ' / '),
      '\n'
    ),
    legText
  )
  legText = paste0(legText,note)
  
  # Samples directory
  version = neutralsEditOrigVersion()
  if(!is.null(neutralsEditCopyVersion()))
    version = neutralsEditCopyVersion()
  neutralsSampleDir = paste0(neutralsPublic,'_',version)
  req(dir.exists(neutralsSampleDir))
  samplesList = list.files(
    path = neutralsSampleDir,
    pattern = 'run_',
    full.names = TRUE
  )
  
  # Generate curves
  T0 = as.numeric(input$T0Plot)
  tRange = seq(input$tempRangePlot[1],input$tempRangePlot[2],5)
  M0  = 10^as.numeric(input$M0Plot)
  mRange = 10^seq(input$densRangePlot[1],input$densRangePlot[2],0.5)
  
  irun = 0
  krateT = matrix(NA,ncol = length(samplesList),nrow= length(tRange))
  krateM = matrix(NA,ncol = length(samplesList),nrow= length(mRange))
  for (file in samplesList) {
    irun = irun + 1
    
    # Get params and generate tags
    scheme  = read.csv(
      file = file,
      skip = iReac-1,
      nrows = 1,
      header = FALSE,
      sep = ';'
    )
    scheme  = gsub(" ", "", scheme)
    params = scheme[8:23]
    pars   = as.numeric(params)
    type   = scheme[ 24]
    
    krateT[,irun] = switch(
      type,
      kooij    = kooij(pars,tRange,M0),
      assocMD  = k3body(pars,tRange,M0),
      assocVV  = kEq18VV(pars,tRange,M0),
      rep(0,length(tRange))
    )
    
    # fixed T, M varies
    krateM[,irun] = switch(
      type,
      kooij    = kooij(pars,T0,mRange),
      assocMD  = k3body(pars,T0,mRange),
      assocVV  = kEq18VV(pars,T0,mRange),
      rep(0,length(mRange))
    )
      
    # if (sum(k <= 0) != 0)
    #   alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
  }
  
  # Plot 
  
  trBlue = col2tr('blue', 60)
  
  par(
    mfrow = c(1, 2),
    mar = c(3,3,9,2),
    mgp = gPars$mgp,
    tcl = gPars$tcl,
    lwd = gPars$lwd,
    pty = 's',
    cex = 1.5
  )
  
  tempRange = tRange
  if (diff(range(krateT)) != 0) {
    matplot(
      tempRange,
      krateT,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr2[6],
      lwd = 1.5*gPars$lwd,
      log = 'y',
      xlab = 'T [K]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    grid()
    lines(tempRange, krateT[, 1], col = gPars$cols[2],lwd = 1.5*gPars$lwd)
    legend('top',title = paste0('M = ',M0,'cm^-3'), legend = NA, bty='n')
    mtext(
      legText,
      side = 3,
      cex = 1.5,
      adj = 0,
      line = 8,
      padj = 1,
      col = gPars$cols[1]
    )
    box()
  }
  
  # P-dep
  if (diff(range(krateM)) != 0) {
    matplot(
      mRange,
      krateM,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr2[6],
      lwd = 1.5*gPars$lwd,
      log = 'xy',
      xlab = 'M [cm^-3]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    lines(mRange, krateM[, 1], col = gPars$cols[2], lwd = 1.5*gPars$lwd)
    grid(col = 'darkgray')
    legend('top',title = paste0('T = ',T0,' K'), legend = NA, bty='n')
    box(lwd = 4)
  }  

},
height = plotHeight)
