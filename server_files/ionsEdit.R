ionsOrigVersion = shiny::reactiveVal()
ionsCopyVersion = shiny::reactiveVal()
ionsFile        = shiny::reactiveVal()
ionsReacs       = shiny::reactiveVal()
ionsRateMask    = shiny::reactiveVal()
ionsBRMask      = shiny::reactiveVal()
ionsSimulSamples= shiny::reactiveVal()

# Generated Inputs ####
output$selIonsOrigVersion = shiny::renderUI({
  list(
    shiny::selectInput(
      "ionsOrigVersion",
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

output$selIonsFile = shiny::renderUI({
  req(ionsOrigVersion())
  list(
    shiny::selectInput(
      "ionsFile",
      "Ions DB File:",
      c("Choose a file..." = "",
        list.files(
          path=file.path(ionsSource,ionsOrigVersion(),'Data'), 
          full.names = FALSE, 
          recursive = FALSE)
      )
    )
  )
})

shiny::observe({
  ionsOrigVersion(input$ionsOrigVersion)
})

shiny::observe({
  ionsCopyVersion(input$ionsCopyVersion)
})

shiny::observe({
  req(input$ionsFile)
  ionsFile(file.path('Data',input$ionsFile,'data.csv'))
})

# Parse data file ####
shiny::observe({
  req(ionsOrigVersion())
  req(ionsFile())
  
  # Get data for text editor
  ionsReacs(
    readLines(
      con <- file(file.path(
        ionsSource,
        ionsOrigVersion(),
        ionsFile()
      ))
    )
  )
  close(con)
  
  # (Re-)Init samples for plots
  ionsSimulSamples(NULL)
  
  # Format data for masks
  mask = list()
  reactants = getSpecies(input$ionsFile)
  mask[['REACTANTS']] = reactants
  massReactants = getMassList(reactants, excludeList = dummySpecies,
                              stoechFilters = stoechFilters)
  # if(length(reactants) > maxReacts)
  #   id = shiny::showNotification(
  #     h4(paste0('Nb of products exceds maxProd=',maxProds)),
  #     closeButton = TRUE,
  #     duration = NULL,
  #     type = 'error'
  #   )
  
  # Get data as matrix
  X = as.matrix(read.csv(
    file.path(ionsSource,ionsOrigVersion(),ionsFile()),
    header=FALSE, sep='\t', fill=TRUE, na.strings=""
  ))
  
  X = trimws(X) # remove unwanted spaces
  reacName = X[1,1]
  if(reacName != input$ionsFile) 
    id = shiny::showNotification(
      h4(paste0('Pb reac identity: ',reacName)),
      closeButton = TRUE,
      duration = NULL,
      type = 'error'
    )
  
  # lastMod = file.info(paste0(dataDir,fName))$mtime
  
  # Locate Rate Info in dataFrame by keywords #
  
  ## Reaction rate expression
  
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
  mask[['TYPE']] = reacType
  
  ## Rate parameters 
  rateParKwdList = c('ALPHA','BETA','GAMMA') 
  for (kwd in rateParKwdList)
    mask[[kwd]] = getDistString(X,kwd)
  
  ## Biblio
  bibKwd = paste0('REF_',c(rateParKwdList,'BR'))
  for (kwd in bibKwd)
    mask[[kwd]] = paste0(getParams(X,kwd),collapse = ';')

  ## Comments
  comments = getParams(X,'RQ')
  if(!is.na(comments)) 
    comments = paste0(comments,collapse=';')
  mask[['RQ']] = comments
  
  ionsRateMask(mask)
  
  # Locate BR Info
  topLeft =  which(X == 'BR', arr.ind = TRUE)
  XBR = X[topLeft[1]:nrow(X), topLeft[2]:ncol(X)]
  
  tags = XBR[2:nrow(XBR), 1] # List of channels
  for (ip in 1:length(tags)) 
    checkBalance(
      getSpecies(input$ionsFile), 
      getSpecies(tags[ip]), 
      stoechFilters = stoechFilters
    )
  ionsBRMask(
    list(
      XBR    = XBR,
      REF_BR = mask[['REF_BR']]
    )
  )
})

output$ionsHeader = shiny::renderUI({
  list(
    h4(file.path(ionsOrigVersion(),ionsFile()))
  )
})

shiny::observe({
  shinyAce::updateAceEditor(
    session, "aceIons", 
    value = paste(ionsReacs(),collapse = '\n')
  )
})

output$ionsRateMask = shiny::renderUI({
  req(ionsRateMask())
  
  mask = ionsRateMask()
  list(
    br(),
    fluidRow(
      column(
        6,
        shiny::textInput(
          "ionsReacReactants",
          "Reactants:",
          value = paste0(mask[['REACTANTS']],collapse = ' + ')
        ),
        shiny::textInput(
          "ionsReacALPHA",
          "ALPHA dist.:",
          value = mask[['ALPHA']]
        ),
        shiny::textInput(
          "ionsReacBETA",
          "BETA dist.",
          value = mask[['BETA']]
        ),
        shiny::textInput(
          "ionsReacGAMMA",
          "GAMMA dist.",
          value = mask[['GAMMA']]
        ),
        shiny::textInput(
          "ionsReacRQ",
          "Comments",
          value = mask[['RQ']]
        )
      ),
      column(
        6,
        shiny::selectInput(
          "ionsReacTYPE",
          "Reaction type:",
          ionsReacTypes,
          selected = mask[['TYPE']]
        ),
        shiny::textInput(
          "ionsReacREF_ALPHA",
          "Refs ALPHA",
          value = mask[['REF_ALPHA']]
        ),
        shiny::textInput(
          "ionsReacREF_BETA",
          "Refs BETA",
          value = mask[['REF_BETA']]
        ),
        shiny::textInput(
          "ionsReacREF_GAMMA",
          "Refs GAMMA",
          value = mask[['REF_GAMMA']]
        )
      )
    )
  )
})
output$ionsBRMask = shiny::renderUI({
  req(ionsRateMask())
  req(ionsBRMask())
  list(
    DT::DTOutput('ionsBR'),
    shiny::textInput(
      "ionsReacREF_BR",
      "Refs BR",
      width = '300px',
      value = ionsBRMask()[['REF_BR']]
    )
  )
})

output$ionsBR = DT::renderDT({
  req(ionsBRMask())
  return(ionsBRMask()[['XBR']])
},
rownames = FALSE,
editable = "cell",
extensions = c('Scroller'),
options = list(
  dom         = 'Btip',
  deferRender = TRUE,
  scrollY     = '500px',
  scroller    = TRUE,
  autoWidth   = TRUE
))

observeEvent(
  input$ionsBR_cell_edit,
  {
    # Apply modifications
    info = input$ionsBR_cell_edit
    XBR = ionsBRMask()$XBR
    for(i in 1:NROW(info))
      XBR[info$row,1+info$col] = info$value
    ionsBRMask(list(XBR=XBR))
  }
)

# Sampling ####
observeEvent(
  input$ionsSimulateBtn,
  {
    req(ionsRateMask)
    req(ionsBRMask)
    
    nMC = as.numeric(input$ionsSampleSize)
    tagged = FALSE 
    
    # Rate parameters
    sampleRateParams = matrix(
      NA, 
      nrow = nMC,
      ncol = length(ionsRateParKwdList))
    colnames(sampleRateParams) = ionsRateParKwdList
    rateParDistStrings = rep(NA,length(ionsRateParKwdList))
    names(rateParDistStrings) = ionsRateParKwdList
    for (kwd in ionsRateParKwdList) {
      stringDist = input[[paste0('ionsReac', kwd)]]
      rateParDistStrings[kwd]=stringDist
      sampleDist = sampleDistString(stringDist, nMC)
      sampleRateParams[1:nMC, kwd] = sampleDist
    }
    
    # Branching ratios
    XBR  = ionsBRMask()$XBR
    tags = XBR[-1,1]
    treeDepth = 0
    for(i in 2:ncol(XBR)) # Exclude NA columns
      treeDepth = treeDepth + 1 - 
      as.numeric(sum(is.na(XBR[,i])) == nrow(XBR))
    
    if (length(tags) >= 2) {
      dist = XBR[1, 2:ncol(XBR)] # List of distributions in tree
      dist = dist[!is.na(dist)]
      XT = XBR[-1, -1]          # Matrix of parameters
      nc = 1
      nl = nrow(XT)
      X1 = matrix(XT[1:nl, 1:nc], nrow = nl, ncol = nc)
      for (ic in 2:ncol(XT))
        if (sum(is.na(XT[, ic])) != nl) {
          nc = nc + 1
          X1 = cbind(X1, XT[, ic])
        }
      
      d = list()
      dSub = list()
      for (ic in 1:nc) {
        dSub$elem = which(X1[, ic] != '')
        dSub$dist = dist[ic]
        pars = as.vector(X1[dSub$elem, ic])
        dSub$mu = dSub$sig = c()
        for (ip in 1:length(pars)) {
          if (grepl("/", pars[ip])) {
            loc  = gregexpr("(?<mu>.*)/(?<sig>.*)", pars, perl = TRUE)
            start = attr(loc[[ip]], 'capture.start')[1]
            stop  = start + attr(loc[[ip]], 'capture.length')[1] - 1
            dSub$mu[ip] = substr(pars[ip], start, stop)
            start = attr(loc[[ip]], 'capture.start')[2]
            stop  = start + attr(loc[[ip]], 'capture.length')[2] - 1
            dSub$sig[ip] = substr(pars[ip], start, stop)
          } else {
            dSub$mu[ip]  = pars[ip]
            dSub$sig[ip] = NA
          }
        }
        dSub$link = rep(0, length(dSub$elem))
        if (ic != 1) {
          for (iel in 1:length(dSub$elem)) {
            ip = dSub$elem[iel]
            for (ic1 in 1:(ic - 1)) {
              links = ip %in% d[[ic1]]$elem
              if (sum(links) != 0)
                dSub$link[iel] = ic1
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
      sampleBR = nds(nMC,stringBR)
      # updateTextInput("ionsStringDist",value=stringBR,session=session)
      
      # Build Newick string for tree plotting 
      newickBR = oneNewick(nc, d, tags)
      while (grepl("LINK/", newickBR)) {
        poc = regmatches(newickBR, gregexpr('LINK/[0-9]+/', newickBR))
        po = sapply(poc[[1]],
                    function(x)
                      as.numeric(sub('LINK', '', gsub('/', '', x))))
        for (ip in 1:length(po)) {
          str = oneNewick(po[ip], d, tags)
          newickBR = sub(poc[[1]][ip], str, newickBR)
        }
      }
      newickBR = paste0(newickBR, ";")
      # updateTextInput(
      #   "ionsStringDist",
      #   value=paste0(newickBR,' ~ ',stringBR),
      #   session=session)
      mytree <- ape::read.tree(text = newickBR)
      
      # Build edge tags for tree annotation #
      edgeTags = oneEdgeTag(nc, d, tags)
      while (grepl("LINK/", edgeTags)) {
        poc = regmatches(edgeTags, gregexpr('LINK/[0-9]+/', edgeTags))
        po = sapply(poc[[1]],
                    function(x)
                      as.numeric(sub('LINK', '', gsub('/', '', x))))
        for (ip in 1:length(po)) {
          str = oneEdgeTag(po[ip], d, tags)
          edgeTags = sub(poc[[1]][ip], str, edgeTags)
        }
      }
      edgeTags = unlist(strsplit(edgeTags, ','))
      
      # Build node tags for tree annotation #
      nodeTags = oneNodeTag(nc, d, tags)
      while (grepl("LINK/", nodeTags)) {
        poc = regmatches(nodeTags, gregexpr('LINK/[0-9]+/', nodeTags))
        po = sapply(poc[[1]],
                    function(x)
                      as.numeric(sub('LINK', '', gsub('/', '', x))))
        for (ip in 1:length(po)) {
          str = oneNodeTag(po[ip], d, tags)
          nodeTags = sub(poc[[1]][ip], str, nodeTags)
        }
      }
      nodeTags = unlist(strsplit(nodeTags, ','))
      nodeTags = nodeTags[nodeTags != 'NA']
      
    } else {
      # Single pathway with BR=1
      sampleBR = matrix(1,ncol=1,nrow=nMC)
      nodeTags = NULL
      edgeTags = NULL
      mytree   = NULL
      treeDepth= NULL
      
    }
    
    ionsSimulSamples(
      list(
        sampleSize       = nMC,
        sampleRateParams = sampleRateParams,
        rateParDistStrings = rateParDistStrings,
        sampleBR         = sampleBR,
        nodeTags         = nodeTags,
        edgeTags         = edgeTags,
        mytree           = mytree,
        treeDepth        = treeDepth,
        tags             = tags
      )
    )
    
  }
)
# Plot rates ####
output$plotIonsParsSample = renderPlot({
  req(ionsSimulSamples())
  req(ionsRateMask())
  
  mask = isolate(ionsRateMask())
  reacType  = mask$TYPE
  reactants = mask$REACTANTS
  
  sampleSize         = ionsSimulSamples()$sampleSize
  sampleRateParams   = ionsSimulSamples()$sampleRateParams
  rateParDistStrings = ionsSimulSamples()$rateParDistStrings
  
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
  
  split.screen(c(2, 1))
  split.screen(c(1,3),screen = 1)
  split.screen(c(1,2),screen = 2)
  iscreen=2
  for (kwd in ionsRateParKwdList) {
    iscreen = iscreen + 1
    screen(iscreen)
    if (!is.finite(sigPars[kwd]) | sigPars[kwd] == 1)
      next
    par(mar = c(4, 5, 3, 1),
        mgp = gPars$mgp,
        tcl = gPars$tcl,
        lwd = gPars$lwd,
        pty = 's',
        cex = 1)
    d = density(sampleRateParams[, kwd])
    plot(
      d,
      xlim = range(sampleRateParams[, kwd]),
      xaxs = 'i',
      xlab = kwd, #paste0(kwd, '~', rateParDistStrings[kwd]),
      yaxs = 'i',
      ylim = c(0,1.1*max(d$y)),
      col = gPars$cols_tr2[5],
      main = ''
    )
    polygon(d$x, d$y, col = gPars$cols_tr[5], border = NA)
    # grid()
    abline(v = meanPars[kwd],
           col = gPars$cols[2],
           lwd = 2)
    abline(
      v = ErrViewLib::vhd(sampleRateParams[, kwd],p = c(0.05,0.95)),
      col = gPars$cols[2],
      lty = 2,
      lwd = 1
    )
    box()
  }
  
  temp = seq(input$ionsTempRangePlot[1],
             input$ionsTempRangePlot[2],
             10)
  nt = length(temp)
  if (reacType == 'kooij') {
    if ('E' %in% reactants) {
      rateFun = function(t, pars)
        pars['ALPHA'] * (300 / t) ^ pars['BETA']
    } else {
      rateFun = function(t, pars)
        pars['ALPHA'] * (t / 300) ^ pars['BETA'] * exp(-pars['GAMMA'] / t)
    }
  } else {
    if (reacType == 'ionpol1') {
      # KIDA::ionpol1 type
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t, pars)
        pars['ALPHA'] * (0.62 + 0.4767 * pars['BETA'] * (300 / t) ^ 0.5)
    } else {
      # KIDA::ionpol2 type
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t, pars)
        pars['ALPHA'] * (1 + 0.0967 * pars['BETA'] * (300 / t) ^ 0.5 +
                           pars['BETA'] ^ 2 * 300 / (10.526 * t))
    }
  }
  
  np = min(sampleSize, 500) # nb of plotted samples
  Y = matrix(NA, nrow = np, ncol = nt)
  for (ip in 1:np)
    Y[ip, 1:nt] = rateFun(temp, sampleRateParams[ip, ])
  yMean = 10^apply(log10(Y),2,mean)
  yF    = 10^apply(log10(Y),2,sd)
  screen(6)
  par(mar = c(4, 5, 1, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1)
  
  matplot(
    temp, t(Y),
    type = 'l',
    lty = 1,
    col = gPars$cols_tr[5],
    lwd = 4,
    log = 'y',
    xlab = 'T / K',
    xaxs = 'i',
    ylab = 'Rate constant / cm^3.s^-1',
    ylim = c(min(Y) / 1.5, max(Y) * 1.5),
    main = ''
  )  
  grid()
  lines(temp,
        yMean, #rateFun(temp, meanPars),
        col = gPars$cols[2],
        lwd = 2)
  box()
  
  screen(7)
  par(mar = c(4, 5, 1, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1)
  plot(
    temp, yF,
    type = 'l',
    lty = 1,
    col = gPars$cols[2],
    lwd = 3,
    xlab = 'T / K',
    xaxs = 'i',
    ylab = 'Uncert. factor',
    main = ''
  )  
  grid()
  box()
  close.screen(all = TRUE)
},
height = plotHeight, width = 1.2*plotWidth)
# Plot BRs ####
## Sample ####
output$plotIonsBRSample = renderPlot({
  req(ionsSimulSamples())
  
  sampleSize       = ionsSimulSamples()$sampleSize
  sampleBR         = ionsSimulSamples()$sampleBR
  tags             = ionsSimulSamples()$tags
  
  meanBR    = colMeans(sampleBR)
  meanBR    = meanBR / sum(meanBR)
  sigBR     = apply(sampleBR, 2, sd)
  
  
  np = min(sampleSize, 500) # nb of plotted samples
  nt = length(tags)
  if (nt >= 2) {
    m0= max(1,40-3*nt)
    par(
      mar = c(m0, 15, 4, 1), 
      cex.main = 1)
    matplot(
      t(sampleBR)[nt:1, 1:np],
      1:nt,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr[5],
      lwd = 5,
      yaxt = 'n',
      ylab = '',
      yaxs = 'i',
      xlim = c(0, 1),
      xlab = '',
      xaxt = 'n',
      xaxs = 'i',
      main = 'Branching Ratios'
    )
    lines(meanBR, nt:1, lwd = 3, col = gPars$cols[2])
    lines(cumsum(meanBR), nt:1, lwd = 3, lty = 2, col = gPars$cols[2])
    text(
      x = seq(0, 1, by = 0.1),
      y = nt,
      labels = seq(0, 1, by = 0.1),
      xpd = TRUE,
      cex = 1,
      pos = 3
    )
    text(
      x = -0.02,
      y = nt:1,
      labels = tags,
      srt = 0,
      adj = 1,
      xpd = TRUE,
      cex = 1
    )
    # grid()
    abline(
      h = 1:length(tags),
      col = 'darkgray',
      lwd = 1,
      lty = 3
    )
    abline(
      v = seq(0.1,1,by=0.1),
      col = 'darkgray',
      lwd = 1,
      lty = 3
    )
    legend('topright', bty = 'n',
           legend = c('Sample', 'Mean', 'Cumul.'),
           lty = c(1,1,2), 
           lwd = 3, pch = NA,
           col = gPars$cols[c(5,2,2)])
    box()
  }
},
height = plotHeight, width = plotWidth)
## Tree ####
output$plotIonsBRTree = renderPlot({
  req(ionsSimulSamples())
  
  sampleSize       = ionsSimulSamples()$sampleSize
  sampleBR         = ionsSimulSamples()$sampleBR
  nodeTags         = ionsSimulSamples()$nodeTags
  edgeTags         = ionsSimulSamples()$edgeTags
  mytree           = ionsSimulSamples()$mytree
  tags             = ionsSimulSamples()$tags
  depth            = ionsSimulSamples()$treeDepth
  
  meanBR    = colMeans(sampleBR)
  meanBR    = meanBR / sum(meanBR)
  sigBR     = apply(sampleBR, 2, sd)
  
  np = min(sampleSize, 500) # nb of plotted samples
  # Branching ratios
  nt = length(tags)
  if (nt >= 2) {
    tagStat = tags
    for (ip in 1:nt) {
      tagStat[ip] = paste0(tags[ip], ' (',
                           signif(meanBR[ip], 2), ' +/- ',
                           signif(sigBR[ip], 1), ')')
    }
    m0 = max(1,40-3*nt)
    m1 = max(1,40-3*depth)
    par(
      mar = c(m0, 0, 1, m1), 
      cex = 1)
    mytree$tip.label = tagStat
    mytree$edge.length = rep(1, dim(mytree$edge)[1])
    plot(
      mytree,
      type = 'clado',
      y.lim = c(length(tags), 1),
      show.tip.label = TRUE,
      tip.color = gPars$cols[1],
      use.edge.length = TRUE,
      root.edge = TRUE,
      edge.width = 2,
      edge.color = gPars$cols[3],
      font = 1,
      main = ''
    )
    mynodelabels(nodeTags, bg = 'gold')
    myedgelabels(edgeTags)
    
  }
},
height = plotHeight, width = 2*plotWidth)

# Biblio ####
output$ionsBiblio = shiny::renderUI({
  req(ionsRateMask())
  req(ionsBRMask())
  
  keys = c()
  for (elt in c('REF_ALPHA','REF_BETA','REF_GAMMA','REF_BR')) {
    k = ionsRateMask()[[elt]]
    if(k != "NA")
      keys = c(keys,unlist(str_split(k,';')))
  }
  keys = sort(unique(keys))
  
  refs = '<H4>References</H4><DL>'
  if(length(keys) != 0) {
    for (key in keys)
      refs = paste0(
        refs,
        '<DT>[',key,']</DT><DD>',
        format(bib[key],style="html"),
        '</DD>')
  }
  refs = paste0(refs,'</DL>')
  list(
    HTML(refs)
  )
})

# output$checkSpecies <- renderText({ 
#   req(input$ace_cursor)
#   sp = input$ace_selection
#   compo = get.atoms(sp, stoechFilters = stoechFilters)
#   names(compo) = elements
#   mass  = massFormula(compo)
#   paste0(
#     "Selection: \"", sp, "\"\n",
#     "Mass= ", mass
#   )
# })

# shiny::observeEvent(
#   input$neuSave,
#   {
#     req(neutralsCopyVersion())
#     req(neutralsFile())
#     
#     # Make copy of source directory
#     neutralsOrigDir = file.path(neutralsSource,neutralsOrigVersion())
#     neutralsCopyDir = file.path(neutralsSource,neutralsCopyVersion())
#     if(!dir.exists(neutralsCopyDir)) {
#       dir.create(neutralsCopyDir)
#       files = list.files(path = neutralsOrigDir)
#       for(file in files)
#         file.copy(
#           from = file.path(neutralsOrigDir,file),
#           to   = file.path(neutralsCopyDir,file)
#         )
#       id = shiny::showNotification(
#         h4(paste0('Created new version: ', neutralsCopyVersion())),
#         closeButton = FALSE,
#         duration = 5
#       )
#     }
#     
#     # Save modified file to target version
#     data = isolate(input$ace)
#     writeLines(
#       data,
#       con = file.path(neutralsCopyDir,neutralsFile())
#     )
#     id = shiny::showNotification(
#       h4(paste0('Saved file: ', neutralsFile())),
#       closeButton = FALSE,
#       duration = 5
#     )
#     
#   }
# )
