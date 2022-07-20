observeEvent(
  input$neutralsSampleBtn,
  {
    req(reacScheme())
    
    # Output directory
    neutralsSampleDir = paste0(neutralsPublic,'_',neutralsVersion())
    if(!dir.exists(neutralsSampleDir))
      dir.create(neutralsSampleDir)
      
    
    nbReac    = reacScheme()$nbReac
    reactants = reacScheme()$reactants
    products  = reacScheme()$products
    params    = reacScheme()$params
    type      = reacScheme()$type
    
    nMC = as.numeric(input$sampleSize)
    
    shiny::withProgress(
      message = 'Sampling...', 
      {
        for (i in 0:nMC) {
          incProgress(1/nMC, detail = if(i==0) 'Nominal' else paste(i,'/',nMC))
          
          if (i==0) # Nominal database
            rnd = rep(0,3*nbReac)
          else
            rnd = rnorm(3*nbReac,0,1)
          
          dbOut = data.frame(
            R1 = 'x', R2 = 'x', R3 = 'x',
            P1 = 'x', P2 = 'x', P3 = 'x', P4 = 'x',
            a = NA,
            b = NA,
            c = NA,
            f = NA,
            g = NA,
            a1 = NA,
            b1 = NA,
            c1 = NA,
            f1 = NA,
            g1 = NA,
            a2 = NA,
            b2 = NA,
            c2 = NA,
            f2 = NA,
            g2 = NA,
            fc = NA,
            ttype = 'x'
          )
          topow = function(x,p) {
            if(x==0)
              return(0.0)
            else
              return(x^p)
          }
          for (m in 1:nbReac) {
            reac = unlist(reactants[m])
            prod = unlist(products[m])
            typ  = type[m]
            pars = unlist(params[m])
            
            if (typ == 'kooij') {
              db1 = data.frame(
                R1 = reac[1], R2 = reac[2], R3 = reac[3],
                P1 = prod[1], P2 = prod[2], P3 = prod[3], P4 = prod[4],
                a = as.numeric(pars[1]),
                b = as.numeric(pars[2]),
                c = as.numeric(pars[3]),
                f = topow(as.numeric(pars[4]),rnd[m]),
                g = as.numeric(pars[5]) * rnd[m],
                a1 = NA,
                b1 = NA,
                c1 = NA,
                f1 = NA,
                g1 = NA,
                a2 = NA,
                b2 = NA,
                c2 = NA,
                f2 = NA,
                g2 = NA,
                fc = NA,
                ttype = typ
              )
            } else {
              if (typ == 'assocMD' | 
                  typ == 'assocVV'  ) {
                db1 = data.frame(
                  R1 = reac[1], R2 = reac[2], R3 = reac[3],
                  P1 = prod[1], P2 = prod[2], P3 = prod[3], P4 = prod[4],
                  a = as.numeric(pars[1]),
                  b = as.numeric(pars[2]),
                  c = as.numeric(pars[3]),
                  f = topow(as.numeric(pars[4]),rnd[m]),
                  g = as.numeric(pars[5]) * rnd[m],
                  a1 = as.numeric(pars[6]),
                  b1 = as.numeric(pars[7]),
                  c1 = as.numeric(pars[8]),
                  f1 = topow(as.numeric(pars[9]),rnd[nbReac + m]),
                  g1 = as.numeric(pars[10]) * rnd[nbReac + m],
                  a2 = as.numeric(pars[11]),
                  b2 = as.numeric(pars[12]),
                  c2 = as.numeric(pars[13]),
                  f2 = topow(as.numeric(pars[14]),rnd[2 * nbReac + m]),
                  g2 = as.numeric(pars[15]) * rnd[2 * nbReac + m],
                  fc = as.numeric(pars[16]),
                  ttype = typ
                )
              }
            }      
            dbOut[m,]=db1
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

output$plotRate = renderPlot({
  req(input$reacNbPlot)
  iReac = as.numeric(input$reacNbPlot)

  req(reacScheme())
  reac0   = reacScheme()$reactants[[iReac]]
  prod0   = reacScheme()$products[[iReac]]
  params0 = reacScheme()$params[[iReac]]
  typ0    = reacScheme()$type[[iReac]]
  note    = reacScheme()$notes[[iReac]]

  tag  = paste0(
    'Reac. ',iReac, ': ',
    paste0(reac0, collapse = ' + '),' --> ',
    paste0(prod0, collapse = ' + ')
  )
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
  neutralsSampleDir = paste0(neutralsPublic,'_',neutralsVersion())
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
    mar = c(3,3,8,2),
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
      line = 7,
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
