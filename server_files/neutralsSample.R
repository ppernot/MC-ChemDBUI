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

  # req(reacScheme())
  # nbReac    = reacScheme()$nbReac
  # reactants = reacScheme()$reactants
  # products  = reacScheme()$products
  # params    = reacScheme()$params
  # type      = reacScheme()$type
  
  # Samples directory
  neutralsSampleDir = paste0(neutralsPublic,'_',neutralsVersion())
  req(dir.exists(neutralsSampleDir))
  samplesList = list.files(
    path = neutralsSampleDir,
    pattern = 'run_',
    full.names = TRUE
  )
  
  # Generate curves
  # T varies, fixed density (M)
  T0 = 150
  tRange = seq(50, 350, by = 5)
  M0  = 1e18 # molec/cm^3
  mRange = 10^seq(8,20,by=1) # molec/cm^3
  
  irun = 0
  krateT = matrix(NA,ncol = length(samplesList),nrow= length(tRange))
  krateM = matrix(NA,ncol = length(samplesList),nrow= length(mRange))
  for (file in samplesList) {
    irun = irun + 1
    
    # Get params and generate tags
    scheme  = read.csv(file = file,
                       skip = iReac-1,
                       nrows = 1,
                       header = FALSE,
                       sep = ';')
    scheme  = gsub(" ", "", scheme)
    terms = scheme[1:3]
    reactants = terms[!is.na(terms) & terms != ""]
    terms = scheme[4:7]
    products = terms[!is.na(terms) & terms != ""]
    terms = scheme[8:23]
    params = terms[!is.na(terms) & terms != ""]
    type = scheme[ 24]
    
    
    pars = as.numeric(params)
    
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
  tag  = paste0(
    iReac,
    ': ',
    paste0(reactants, collapse =
             '+'),
    '->',
    paste0(products, collapse = '+')
  )
  legText = paste0(tag, '\n',
                   'Rate law: ', type)
  
  legText = switch(
    type,
    kooij    = paste0(
      legText,
      '\n',
      'Parameters: ',
      paste0(params[1:5],
             collapse = ' / '),
      '\n'
    ),
    assocMD  = paste0(
      legText,
      '\n',
      'Parameters Fc : ',
      params[16],
      '\n',
      'Parameters k0 : ',
      paste0(params[1:5], collapse = ' / '),
      '\n',
      'Parameters kInf : ',
      paste0(params[6:10], collapse = ' / '),
      '\n',
      'Parameters kr : ',
      paste0(params[11:15], collapse = ' / '),
      '\n'
    ),
    assocVV  = paste0(
      legText,
      '\n',
      'Parameters Fc : ',
      params[16],
      '\n',
      'Parameters kInf : ',
      paste0(params[1:5], collapse = ' / '),
      '\n',
      'Parameters k0 : ',
      paste0(params[6:10], collapse = ' / '),
      '\n',
      'Parameters kR : ',
      paste0(params[11:15], collapse = ' / '),
      '\n'
    ),
    legText
  )
  
  col2tr = function(x, alpha = 80) {
    rgb(unlist(t(col2rgb(x))), alpha = alpha, maxColorValue = 255)
  }
  trBlue = col2tr('blue', 60)
  
  par(
    mfrow = c(1, 2),
    mar = c(4, 4, 10, 2),
    cex = 1
  )
  
  tempRange = tRange
  if (diff(range(krateT)) != 0) {
    matplot(
      tempRange,
      krateT,
      type = 'l',
      lty = 1,
      col = trBlue,
      lwd = 3,
      log = 'y',
      xlab = 'T [K]',
      ylab = 'rate ct. [cm^3.s^-1]',
      main = ''
    )
    lines(tempRange, krateT[, 1], col = 'red', lwd = 3)
    grid(col = 'darkgray')
    legend('topright',title = 'M = 1e18 cm^-3', legend = NA, bty='n')
    mtext(
      legText,
      side = 3,
      cex = 1.5,
      adj = 0,
      line = 9,
      padj = 1,
      col = 'darkgreen'
    )
    box(lwd = 4)
  }
  
  # P-dep
  if (diff(range(krateM)) != 0) {
    matplot(
      mRange,
      krateM,
      type = 'l',
      lty = 1,
      col = trBlue,
      lwd = 3,
      log = 'xy',
      xlab = 'M [cm^-3]',
      ylab = 'rate ct. [cm^3.s^-1]',
      main = ''
    )
    lines(mRange, krateM[, 1], col = 'red', lwd = 3)
    grid(col = 'darkgray')
    legend('topright',title = 'T = 150 K', legend = NA, bty='n')
    box(lwd = 4)
  }  

},
height = plotHeight)
