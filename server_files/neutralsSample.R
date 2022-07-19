observeEvent(
  input$neutralsSampleBtn,
  {
    req(reacScheme())
    
    nbReac    = reacScheme()$nbReac
    reactants = reacScheme()$reactants
    products  = reacScheme()$products
    params    = reacScheme()$params
    type      = reacScheme()$type
    
    withProgress(
      message = 'Reactions', 
      {
        
        for (i in 0:input$sampleSize) {
          
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
          # write.table(dbOut,
          #             file=paste0(samplesDir,'run_',sprintf('%04i',i),'.csv'),
          #             quote=TRUE, sep=';', na=" ",
          #             row.names=FALSE,col.names=FALSE)  
          
          incProgress(1/input$sampleSize, detail = paste('Sample', i))
        }
        
      })
    
  })

