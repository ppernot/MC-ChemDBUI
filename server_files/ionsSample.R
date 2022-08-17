collateDone = shiny::reactiveVal(NULL)

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
# Sampling ####
observeEvent(
  input$ionsSampleBtn,
  {
    
    sampleSize = 1 + as.numeric(input$ionsSampleSize) # Account for nominal
    tagged = FALSE 
    
    # List of reacs in DB
    fp        = file.path(ionsSource,input$ionsVersionSample,'Data')
    listDirs  = list.files(path = fp, full.names = FALSE, recursive = FALSE)
    tmpDir    = file.path(ionsTmp,input$ionsVersionSample)
    fpTmp     = file.path(tmpDir,'Reactions')
    
    ## Sampling loop ####
    shiny::withProgress(
      message = 'Sampling ', 
      {
        sampled = c()
        for (reac in listDirs) {
          
          # Check if sampling is necessary
          needSample = FALSE
          if(input$ionsSampleUpdate) {
            # Filter by modification date and sample size
            dataFile  = file.path(fp,reac,'data.csv') 
            finf1     = file.info(dataFile, extra_cols = FALSE) # Data mod times
            
            targetDir = file.path(fpTmp,reac)
            if(!file.exists(targetDir)) {
              dir.create(targetDir)
              needSample = TRUE
            }
            
            samplesDir = file.path(targetDir, 'Samples')
            if (!file.exists(samplesDir)){
              dir.create(sampleDir)
              needSample = TRUE
            }
            
            tmpSampleSize = length(list.files(samplesDir))
            if(tmpSampleSize < input$ionsSampleSize)
              needSample = TRUE
            
            runFile = file.path(samplesDir,'run_0000.csv')
            if(file.exists(runFile)) {
              finf2     = file.info(runFile, extra_cols = FALSE) # Samples mod times
              if(difftime(finf1[,"mtime"],finf2[,"mtime"]) > 0 )
                needSample = TRUE
            } else{
              needSample = TRUE
            }
          }
          
          incProgress(1/length(listDirs), detail = reac)
          
          if(!needSample) next()
          
          sampled = c(sampled, reac)
          
          allBibKeys = allSpecies = c()
          
          reactants = getSpecies(reac)
          allSpecies = c(allSpecies,reactants)
          massReactants = getMassList(reactants, excludeList = dummySpecies)
          
          # Get data for this reaction #
          X = as.matrix(
            read.csv(
              file.path(fp,reac,'data.csv'),
              header=FALSE, sep='\t', fill=TRUE, na.strings=""
            )
          )
          
          X = trimws(X) # remove unwanted spaces
          reacName = X[1,1]
          if(reacName != reac) 
            id = shiny::showNotification(
              h4(paste0('Pb reac identity: ',reacName)),
              closeButton = TRUE,
              duration = NULL,
              type = 'error'
            )
          
          # Locate Rate Info in X by keywords #
          topLeft = which(X=='TYPE',arr.ind=TRUE)
          if(length(topLeft)==0) {
            reacType = 'kooij' # default
            if('E' %in% reactants) reacType ='dr'
          } else {  
            reacType = X[topLeft[1],topLeft[2]+1]
            if(! reacType %in% ionsReacTypes) 
              id = shiny::showNotification(
                h4(paste0('Improper rate type:',reacType)),
                closeButton = TRUE,
                duration = NULL,
                type = 'error'
              )
          }
          
          # Rate parameters 
          sampleRateParams = matrix(NA,
                                    nrow = sampleSize,
                                    ncol = length(ionsRateParKwdList))
          colnames(sampleRateParams) = ionsRateParKwdList
          rateParDistStrings = rep(NA, length(ionsRateParKwdList))
          names(rateParDistStrings) = ionsRateParKwdList
          for (kwd in ionsRateParKwdList) {
            stringDist = getDistString(X, kwd)
            rateParDistStrings[kwd] = stringDist
            sampleDist = sampleDistString(stringDist, sampleSize)
            sampleRateParams[1:sampleSize, kwd] = sampleDist
          }
          
          # Locate BR Info #
          topLeft =  which(X == 'BR', arr.ind = TRUE)
          XBR = X[topLeft[1]:nrow(X), topLeft[2]:ncol(X)]
          tags = XBR[2:nrow(XBR), 1] # List of channels
          
          # Check mass compatibility
          allProds = c()
          for (ip in 1:length(tags)) {
            prods = getSpecies(tags[ip])
            allProds = c(allProds, prods)
            allSpecies = c(allSpecies, prods)
            msg = checkBalance(reactants, prods)
            if (!is.null(msg))
              id = shiny::showNotification(h4(msg),
                                           closeButton = TRUE,
                                           duration = NULL,
                                           type = 'error')
          }
          
          if (input$ionsSampleCheck) next
          
          ### Generate BR samples ####
          if(length(tags) >=2) {
            # Build tree
            dist=XBR[1,2:ncol(XBR)] # List of distributions in tree
            dist=dist[!is.na(dist)]
            XT=XBR[-1,-1]          # Matrix of parameters
            nc=1
            nl=nrow(XT)
            X1=matrix(XT[1:nl,1:nc],nrow=nl,ncol=nc)
            for (ic in 2:ncol(XT)) 
              if(sum(is.na(XT[,ic]))!=nl) {
                nc=nc+1
                X1=cbind(X1,XT[,ic])
              }
            
            d=list()
            dSub=list()
            for(ic in 1:nc){
              dSub$elem = which(X1[,ic] !='')
              dSub$dist = dist[ic]
              pars = as.vector(X1[dSub$elem,ic])
              dSub$mu = dSub$sig = c()
              for (ip in 1:length(pars)) {
                if( grepl("/",pars[ip]) ) {      
                  loc  = gregexpr("(?<mu>.*)/(?<sig>.*)",pars,perl=TRUE)      
                  start = attr(loc[[ip]],'capture.start')[1]
                  stop  = start + attr(loc[[ip]],'capture.length')[1] -1
                  dSub$mu[ip] = substr(pars[ip],start,stop)
                  start = attr(loc[[ip]],'capture.start')[2]
                  stop  = start + attr(loc[[ip]],'capture.length')[2] -1
                  dSub$sig[ip] = substr(pars[ip],start,stop)  
                } else {
                  dSub$mu[ip]  = pars[ip]
                  dSub$sig[ip] = NA
                }    
              }
              dSub$link = rep(0,length(dSub$elem))
              if(ic!=1) {
                for (iel in 1:length(dSub$elem)) {
                  ip = dSub$elem[iel]
                  for (ic1 in 1:(ic-1)) {
                    links = ip %in% d[[ic1]]$elem
                    if(sum(links)!=0) dSub$link[iel] = ic1  
                  }
                }
              }
              d[[ic]] = dSub
            }
            
            # Build probabilistic tree string for sampler #
            stringBR = oneDist(nc, d, tags, tagged)
            while (grepl("LINK/", stringBR)) {
              poc = regmatches(stringBR, gregexpr('LINK/[0-9]+/', stringBR))
              po = sapply(poc[[1]],
                          function(x)
                            as.numeric(sub('LINK', '', gsub('/', '', x))))
              for (ip in 1:length(po)) {
                str = oneDist(po[ip], d, tags, tagged)
                stringBR = sub(poc[[1]][ip], str, stringBR)
              }
            }    
            # Generate BR sample #
            sampleBR = nds(sampleSize,stringBR)  
            
          } else {
            # Single pathway with BR=1
            sampleBR = matrix(1,ncol=1,nrow=sampleSize)
            
          }
          
          # Generate output for kinetics code #
          
          # Nominal/mean/median values from samples
          meanPars = rep(NA, ncol(sampleRateParams))
          names(meanPars) = colnames(sampleRateParams)
          sigPars = rep(NA, ncol(sampleRateParams))
          names(sigPars) = colnames(sampleRateParams)
          for (kwd in ionsRateParKwdList) {
            sample = sampleRateParams[, kwd]
            if(substr(rateParDistStrings[kwd],1,5)!='Delta') {
              meanPars[kwd] = exp(mean(log(sample)))
              sigPars[kwd]  = exp(sd(log(sample)))
            } else {
              meanPars[kwd] = mean(sample)
              sigPars[kwd]  = 1
            }
          }
          rm(sample)
          
          meanBR    = colMeans(sampleBR)
          meanBR    = meanBR / sum(meanBR)
          sigBR     = apply(sampleBR, 2, sd)
          
          
          ## Write samples ####
          # Write sample to tmp file
          writeSample(0, samplesDir, reac, tags, meanPars, meanBR, reacType)
          for (i in 1:sampleSize) 
            writeSample(i, samplesDir, reac, tags, 
                        sampleRateParams[i,], sampleBR[i,], 
                        reacType)     
          
          # # Auxiliary files
          # sink(file = file.path(targetDir, 'species.txt'), append = FALSE)
          # cat(unique(allSpecies))
          # sink(file = NULL)
          
          bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
          for (kwd in bibKwd) {
            refs = getParams(X,kwd)
            allBibKeys = c(allBibKeys,refs) 
          }
          allBibKeys=sort(unique(allBibKeys[!is.na(allBibKeys)]))
          
          sink(file = file.path(targetDir, 'bibKeys.txt'), append = FALSE)
          cat(allBibKeys)
          sink(file = NULL)
          
          sink(file=file.path(targetDir,'dataTable.html'),append=FALSE)
          cat('<TR>\n')
          cat(paste0('<TD>',reac,'</TD>\n'))
          cat(paste0('<TD>--></TD>\n'))
          cat(paste0('<TD>',tags[1],'</TD>\n'))
          if(length(tags) <2) {
            mbr = 1
            sbr = 0
          } else{
            mbr = meanBR[1]
            sbr = sigBR[1]    
          }
          cat(paste0('<TD>',paste0(signif(mbr,2),' +/- ',
                                   signif(sbr,1))    ,'</TD>\n'))
          cat(paste0('<TD>',paste0(signif(meanPars['ALPHA'],2),' */ ',
                                   signif(sigPars['ALPHA'],2))    ,'</TD>\n'))
          cat('</TR>\n')
          if(length(tags) >= 2) {
            for(i in 2:length(tags)) {
              cat('<TR>\n')
              cat(paste0('<TD> </TD>\n'))
              cat(paste0('<TD> </TD>\n'))
              cat(paste0('<TD>',tags[i],'</TD>\n'))
              cat(paste0('<TD>',paste0(signif(meanBR[i],2),' +/- ',
                                       signif(sigBR[i],1))    ,'</TD>\n'))
              cat(paste0('<TD> </TD>\n'))
              cat('</TR>\n')      
            }
          }
          sink(file = NULL)
          
        }
      }
    )
    if(length(sampled)== 0) {
      id = shiny::showNotification(
        h4('Sampling: nothing to be done...'),
        closeButton = TRUE,
        duration = NULL,
        type = 'warning'
      )
    } else if(input$ionsSampleUpdate) {
      id = shiny::showNotification(
        h4(paste0('Updated ',paste0(sampled,collapse = ', '))),
        closeButton = TRUE,
        duration = NULL,
        type = 'warning'
      )
    }
    
    ## Gather loop ####
    
    if (!input$ionsSampleCheck) {
      
      # List of reacs in Tmp DB
      tmpDir    = file.path(ionsTmp,input$ionsVersionSample,'Reactions')
      listReacs = list.files(path = tmpDir, full.names = FALSE, recursive = FALSE)
      
      # Target directory
      ionsSampleDir = paste0(ionsPublic, '_', input$ionsVersionSample)
      if (!dir.exists(ionsSampleDir)) {
        dir.create(ionsSampleDir)
      } else {
        # Clean
        fileList = list.files(path = ionsSampleDir, full.name = TRUE)
        file.remove(fileList)
      }
      
      dataTableFile = file.path(ionsSampleDir,'dataTable.html')
      if(file.exists(dataTableFile))
        file.remove(dataTableFile)
      sink(file=dataTableFile)
      cat('<!DOCTYPE html>\n<HTML>\n
      <HEAD>\n<STYLE>th {text-align: left;}</STYLE>\n</HEAD>\n
      <BODY><TABLE BORDER=0>\n
      <TR><TH>Reactants</TH><TH></TH><TH>Products</TH>
          <TH>Branching Ratios</TH><TH>Alpha</TH></TR>\n')
      sink(file=NULL)
      
      shiny::withProgress(
        message = 'Collating ', 
        {
          allSpecies = allBibKeys = c()
          for (reac in listReacs) {
            
            incProgress(1/length(listReacs), detail = reac)
            
            sink(file = dataTableFile, append = TRUE)
            cat('<TR><TD COLSPAN=5><HR size=1></TD></TR>\n')
            sink(file = NULL)
            file.append(
              file1 = dataTableFile,
              file2 = file.path(tmpDir, reac, 'dataTable.html')
            )
            
            # Generate collated Monte Carlo samples
            for (i in 0:(sampleSize-1)) {
              runFile = paste0('run_', sprintf('%04i', i), '.csv')
              file.append(
                file1 = file.path(ionsSampleDir, runFile),
                file2 = file.path(tmpDir, reac, 'Samples', runFile)
              )
            }
            
            # Collate full biblio
            file = file.path(tmpDir, reac, 'bibKeys.txt')
            if (file.info(file)$size != 0) {
              bibLine = readLines(con = file(file), n = 1)
              bibKeys = str_split(bibLine, ' ')
              allBibKeys = c(allBibKeys, unlist(bibKeys))
            }
            
          } 
        })
      
      sink(file=dataTableFile,append=TRUE)
      cat('</TABLE></HTML>')
      sink(file=NULL)
      
      targetHtml = file.path(ionsSampleDir,'bibliography.html')
      sink(file=targetHtml, append=FALSE)
      allBibKeys=sort(unique(allBibKeys))
      printBib(allBibKeys,bib)
      sink(file = NULL)
      
      id = shiny::showNotification(
        h4('Samples written to ChemDBPublic'),
        closeButton = TRUE,
        duration = NULL,
        type = 'message'
      )
      
      collateDone(TRUE)
    }   
  })

# Statistics ####
output$ionsStats = renderPrint({
  req(collateDone())
  
  id = shiny::showNotification(
    h4('Computing statistics...'),
    closeButton = TRUE,
    duration = NULL,
    type = 'message'
  )
  
  fileName = file.path(
    paste0(ionsPublic,'_',input$ionsVersionSample),'run_0000.csv')
  scheme   = as.data.frame(
    data.table::fread(file = fileName, header = FALSE, sep = ';')
  )
  nbReac   = nrow(scheme)
  reactants = products = params = type = reacTagFull = list()
  for (i in 1:nbReac) {
    line = replace(scheme[i, ], scheme[i, ] == " ", "")
    line = replace(line, is.na(line), "")
    reactants[[i]] = line[1:3][line[1:3]!=""]
    products[[i]]  = line[4:7][line[4:7]!=""]
    type[[i]]      = line[ncol(scheme)]
    reacTagFull[[i]] = paste0(
      paste(reactants[[i]], collapse = ' + '), 
      ' --> ', 
      paste(products[[i]], collapse = ' + ')
    )
  }
  
  # Build species list from reactants and products
  species   = unique(unlist(c(reactants,products)))
  nbSpecies = length(species)
  reacts    = unique(unlist(reactants))
  nbReacts  = length(reacts)
  prods     = unique(unlist(products))
  nbProds   = length(prods)
  
  cat('Nb reactions       = ', nbReac   ,'\n')
  cat('Nb species         = ', nbSpecies,'\n\n')
  
  # Nature of species
  types = c("neutrals","ions")
  chems = c("hydrocarbons","N-bearing","O-bearing")
  heavy = c("C0","C1","C2","C3","C4","C5","C6","Cmore")
  resu = matrix(0,
                nrow = 1 + length(chems),
                ncol = 2 + length(types))
  rownames(resu) = c(chems,'total')
  colnames(resu) = c(types,'radicals','total')
  for (type in types)
    resu['total',type] = sum(
      selectSpecies(species,c(type,chems,heavy)))
  resu['total','radicals'] = sum(
    selectSpecies(species,c(types,'radicals',chems,heavy)))
  for (chem in chems) {
    resu[chem,'total'] = sum(
      selectSpecies(species,c(types,chem,heavy)))
    resu[chem,'radicals'] = sum(
      selectSpecies(species,c(types,'radicals',chem,heavy)))
  }
  for (type in types)
    for (chem in chems)
      resu[chem,type] = sum(
        selectSpecies(species,c(type,chem,heavy)))
  resu['total','total'] = sum(resu[,'total'])
  cat('Compositions (excl. dummies)\n')
  cat('----------------------------\n\n')
  print(resu)
  cat('\n\n')
  
  # Loss and prod
  L = R = matrix(0,ncol=nbSpecies,nrow=nbReac)
  for (m in 1:nbReac) {
    reac = unlist(reactants[m])
    prod = unlist(products[m] )
    for (n in 1:nbSpecies) {
      search=species[n]
      L[m,n] = length(which( search == reac )) # Loss
      R[m,n] = length(which( search == prod )) # Prod
    }
  }
  selLossless = colSums(L) == 0
  nbLossless = sum(selLossless)
  if(nbLossless != 0)
    lossless = sort(species[selLossless])
  selProdless = colSums(R) == 0
  nbProdless = sum(selProdless)
  if(nbProdless != 0)
    prodless = sort(species[selProdless])

  cat('Reactants   : ',nbReacts,'\n')
  cat('Products    : ',nbProds,'\n\n')
  cat('Number of lossless species = ',nbLossless)
  if(nbLossless != 0) 
    cat(strwrap(paste0(lossless, collapse = ', '),width=60,prefix='\n'),'\n\n')
  cat('Number of prodless species = ',nbProdless)
  if(nbProdless != 0) 
    cat(strwrap(paste0(prodless, collapse = ', '),width=60,prefix='\n'),'\n\n')
  
  id0=shiny::removeNotification(id)
})


# Convert to single-file DB ####
observeEvent(
  input$ionsCnvrtBtn,
  {
    
    # List of reacs in DB
    fp        = file.path(ionsSource,input$ionsVersionSample,'Data')
    listDirs  = list.files(path = fp, full.names = FALSE, recursive = FALSE)
    
    ionsSampleDir = paste0(ionsPublic, '_', input$ionsVersionSample)
    outFile = file.path(ionsSampleDir,'ionsDB.csv')
    if(file.exists(outFile))
      file.remove(outFile)
    
    bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
    sink(file = outFile, append = TRUE)
    cat(
      'REACTANTS;',
      'TYPE;',
      paste0(ionsRateParKwdList,collapse=';'),';',
      'NBR;',
      'STRINGBR;',
      paste0(bibKwd,collapse=';'),';',
      'COMMENTS\n'
    )
    sink(file = NULL)

    ## Reactions loop ####
    shiny::withProgress(
      message = 'Converting ', 
      {
        allSpecies = c()
        sampled = c()
        for (reac in listDirs) {
          
          incProgress(1/length(listDirs), detail = reac)
          
          sampled = c(sampled, reac)
          
          reactants = getSpecies(reac)
          
          # Get data for this reaction #
          X = as.matrix(
            read.csv(
              file.path(fp,reac,'data.csv'),
              header=FALSE, sep='\t', fill=TRUE, na.strings=""
            )
          )
          
          X = trimws(X) # remove unwanted spaces
          reacName = X[1,1]
          if(reacName != reac) 
            id = shiny::showNotification(
              h4(paste0('Pb reac identity: ',reacName)),
              closeButton = TRUE,
              duration = NULL,
              type = 'error'
            )
          
          # Locate Rate Info in X by keywords #
          topLeft = which(X=='TYPE',arr.ind=TRUE)
          if(length(topLeft)==0) {
            reacType = 'kooij' # default
            if('E' %in% reactants) reacType ='dr'
          } else {  
            reacType = X[topLeft[1],topLeft[2]+1]
            if(! reacType %in% ionsReacTypes) 
              id = shiny::showNotification(
                h4(paste0('Improper rate type:',reacType)),
                closeButton = TRUE,
                duration = NULL,
                type = 'error'
              )
          }
          
          # Rate parameters 
          rateParDistStrings = rep(NA, length(ionsRateParKwdList))
          names(rateParDistStrings) = ionsRateParKwdList
          for (kwd in ionsRateParKwdList) {
            stringDist = getDistString(X, kwd)
            rateParDistStrings[kwd] = stringDist
          }
          
          refBib = rep(NA, length(bibKwd))
          names(refBib) = bibKwd
          for (kwd in bibKwd)
            refBib[kwd] = paste0('"',paste0(getParams(X,kwd),collapse = ';'),'"')
          
          comments = getParams(X,'RQ')
          if(!is.na(comments)) 
            comments = paste0(comments,collapse=';')
          
          # Locate BR Info #
          topLeft =  which(X == 'BR', arr.ind = TRUE)
          XBR = X[topLeft[1]:nrow(X), topLeft[2]:ncol(X)]
          tags = XBR[2:nrow(XBR), 1] # List of channels
          
          # Check mass compatibility
          allProds = c()
          for (ip in 1:length(tags)) {
            prods = getSpecies(tags[ip])
            allProds = c(allProds, prods)
            allSpecies = c(allSpecies, prods)
            msg = checkBalance(reactants, prods)
            if (!is.null(msg))
              id = shiny::showNotification(h4(msg),
                                           closeButton = TRUE,
                                           duration = NULL,
                                           type = 'error')
          }
          
          ### Generate BR distrib ####
          if(length(tags) >=2) {
            # Build tree
            dist=XBR[1,2:ncol(XBR)] # List of distributions in tree
            dist=dist[!is.na(dist)]
            XT=XBR[-1,-1]          # Matrix of parameters
            nc=1
            nl=nrow(XT)
            X1=matrix(XT[1:nl,1:nc],nrow=nl,ncol=nc)
            for (ic in 2:ncol(XT)) 
              if(sum(is.na(XT[,ic]))!=nl) {
                nc=nc+1
                X1=cbind(X1,XT[,ic])
              }
            
            d=list()
            dSub=list()
            for(ic in 1:nc){
              dSub$elem = which(X1[,ic] !='')
              dSub$dist = dist[ic]
              pars = as.vector(X1[dSub$elem,ic])
              dSub$mu = dSub$sig = c()
              for (ip in 1:length(pars)) {
                if( grepl("/",pars[ip]) ) {      
                  loc  = gregexpr("(?<mu>.*)/(?<sig>.*)",pars,perl=TRUE)      
                  start = attr(loc[[ip]],'capture.start')[1]
                  stop  = start + attr(loc[[ip]],'capture.length')[1] -1
                  dSub$mu[ip] = substr(pars[ip],start,stop)
                  start = attr(loc[[ip]],'capture.start')[2]
                  stop  = start + attr(loc[[ip]],'capture.length')[2] -1
                  dSub$sig[ip] = substr(pars[ip],start,stop)  
                } else {
                  dSub$mu[ip]  = pars[ip]
                  dSub$sig[ip] = NA
                }    
              }
              dSub$link = rep(0,length(dSub$elem))
              if(ic!=1) {
                for (iel in 1:length(dSub$elem)) {
                  ip = dSub$elem[iel]
                  for (ic1 in 1:(ic-1)) {
                    links = ip %in% d[[ic1]]$elem
                    if(sum(links)!=0) dSub$link[iel] = ic1  
                  }
                }
              }
              d[[ic]] = dSub
            }
            
            # Build probabilistic tree string for sampler #
            stringBR = oneDist(nc, d, tags, tagged=TRUE)
            while (grepl("LINK/", stringBR)) {
              poc = regmatches(stringBR, gregexpr('LINK/[0-9]+/', stringBR))
              po = sapply(poc[[1]],
                          function(x)
                            as.numeric(sub('LINK', '', gsub('/', '', x))))
              for (ip in 1:length(po)) {
                str = oneDist(po[ip], d, tags, tagged=TRUE)
                stringBR = sub(poc[[1]][ip], str, stringBR)
              }
            }    

          } else {
            stringBR = tags
          }
          
          sink(file = outFile, append = TRUE)
          cat(
            paste0(
              '"',trimws(reac),'";',
              '"',trimws(reacType),'";',
              paste0('"',rateParDistStrings,'"',collapse=';'),';',
              length(tags),';',
              '"',trimws(stringBR),'";',
              paste0(refBib,collapse=';'),';',
              '"',comments,'"\n'
            )
          )
          sink(file = NULL)
          
        }
      }
    )

    id = shiny::showNotification(
      h4('Samples written to ChemDBPublic'),
      closeButton = TRUE,
      duration = NULL,
      type = 'message'
    )
  })
