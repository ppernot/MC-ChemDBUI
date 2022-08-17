# Modif of ape function to enlarge "rect" frames in labels
mynodelabels = function (text,
                         node,
                         adj = c(0.5, 0.5),
                         frame = "rect",
                         pch = NULL,
                         thermo = NULL,
                         pie = NULL,
                         piecol = NULL,
                         col = "black",
                         bg = "lightblue",
                         horiz = FALSE,
                         width = NULL,
                         height = NULL,
                         ...)
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(node))
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
  XX <- lastPP$xx[node]
  YY <- lastPP$yy[node]
  myBOTHlabels(
    text,
    node,
    XX,
    YY,
    adj,
    frame,
    pch,
    thermo,
    pie,
    piecol,
    col,
    bg,
    horiz,
    width,
    height,
    ...
  )
}
myedgelabels = function (text,
                         edge,
                         adj = c(0.5, 0.5),
                         frame = "rect",
                         pch = NULL,
                         thermo = NULL,
                         pie = NULL,
                         piecol = NULL,
                         col = "black",
                         bg = "lightgreen",
                         horiz = FALSE,
                         width = NULL,
                         height = NULL,
                         date = NULL,
                         ...)
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(edge)) {
    sel <- 1:dim(lastPP$edge)[1]
    subedge <- lastPP$edge
  }
  else {
    sel <- edge
    subedge <- lastPP$edge[sel, , drop = FALSE]
  }
  if (lastPP$type == "phylogram") {
    if (lastPP$direction %in% c("rightwards", "leftwards")) {
      XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[,
                                                         2]]) / 2
      YY <- lastPP$yy[subedge[, 2]]
    }
    else {
      XX <- lastPP$xx[subedge[, 2]]
      YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[,
                                                         2]]) / 2
    }
  }
  else {
    XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[,
                                                       2]]) / 2
    YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[,
                                                       2]]) / 2
  }
  if (!is.null(date))
    XX[] <- max(lastPP$xx) - date
  myBOTHlabels(text,
               sel,
               XX,
               YY,
               adj,
               frame,
               pch,
               thermo,
               pie,
               piecol,
               col,
               bg,
               horiz,
               width,
               height,
               ...)
}
myBOTHlabels = function (text,
                         sel,
                         XX,
                         YY,
                         adj,
                         frame,
                         pch,
                         thermo,
                         pie,
                         piecol,
                         col,
                         bg,
                         horiz,
                         width,
                         height,
                         ...)
{
  if (missing(text))
    text <- NULL
  if (length(adj) == 1)
    adj <- c(adj, 0.5)
  if (is.null(text) &&
      is.null(pch) && is.null(thermo) && is.null(pie))
    text <- as.character(sel)
  frame <- match.arg(frame, c("rect", "circle", "none"))
  args <- list(...)
  CEX <- if ("cex" %in% names(args))
    args$cex
  else
    par("cex")
  if (frame != "none" && !is.null(text)) {
    if (frame == "rect") {
      width <- strwidth(text, units = "inches", cex = CEX)
      height <- strheight(text, units = "inches", cex = CEX)
      if ("srt" %in% names(args)) {
        args$srt <- args$srt %% 360
        if (args$srt == 90 || args$srt == 270) {
          tmp <- width
          width <- height
          height <- tmp
        }
        else if (args$srt != 0)
          warning(
            "only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n"
          )
      }
      width <- xinch(width)
      height <- yinch(height)
      xl <- XX - width * adj[1] - xinch(0.04)
      xr <- xl + width + xinch(0.08)
      yb <- YY - height * adj[2] - yinch(0.04)
      yt <- yb + height + yinch(0.08)
      rect(xl, yb, xr, yt, col = bg)
    }
    if (frame == "circle") {
      radii <- 0.8 * apply(cbind(
        strheight(text, units = "inches",
                  cex = CEX),
        strwidth(text, units = "inches",
                 cex = CEX)
      ), 1, max)
      symbols(
        XX,
        YY,
        circles = radii,
        inches = max(radii),
        add = TRUE,
        bg = bg
      )
    }
  }
  if (!is.null(thermo)) {
    parusr <- par("usr")
    if (is.null(width)) {
      width <- CEX * (parusr[2] - parusr[1])
      width <- if (horiz)
        width / 15
      else
        width / 40
    }
    if (is.null(height)) {
      height <- CEX * (parusr[4] - parusr[3])
      height <- if (horiz)
        height / 40
      else
        height / 15
    }
    if (is.vector(thermo))
      thermo <- cbind(thermo, 1 - thermo)
    thermo <- if (horiz)
      width * thermo
    else
      height * thermo
    if (is.null(piecol))
      piecol <- rainbow(ncol(thermo))
    xl <- XX - width / 2 + adj[1] - 0.5
    xr <- xl + width
    yb <- YY - height / 2 + adj[2] - 0.5
    yt <- yb + height
    if (horiz) {
      rect(xl,
           yb,
           xl + thermo[, 1],
           yt,
           border = NA,
           col = piecol[1])
      for (i in 2:ncol(thermo))
        rect(
          xl + rowSums(thermo[,
                              1:(i - 1), drop = FALSE]),
          yb,
          xl + rowSums(thermo[,
                              1:i]),
          yt,
          border = NA,
          col = piecol[i]
        )
    }
    else {
      rect(xl,
           yb,
           xr,
           yb + thermo[, 1],
           border = NA,
           col = piecol[1])
      for (i in 2:ncol(thermo))
        rect(
          xl,
          yb + rowSums(thermo[,
                              1:(i - 1), drop = FALSE]),
          xr,
          yb + rowSums(thermo[,
                              1:i]),
          border = NA,
          col = piecol[i]
        )
    }
    s <- apply(thermo, 1, function(xx)
      any(is.na(xx)))
    xl[s] <- xr[s] <- NA
    rect(xl, yb, xr, yt, border = "black")
    if (!horiz) {
      segments(xl, YY, xl - width / 5, YY)
      segments(xr, YY, xr + width / 5, YY)
    }
  }
  if (!is.null(pie)) {
    if (is.vector(pie))
      pie <- cbind(pie, 1 - pie)
    xrad <- CEX * diff(par("usr")[1:2]) / 50
    xrad <- rep(xrad, length(sel))
    XX <- XX + adj[1] - 0.5
    YY <- YY + adj[2] - 0.5
    for (i in seq_along(sel)) {
      if (any(is.na(pie[i,])))
        next
      floating.pie.asp(XX[i], YY[i], pie[i,], radius = xrad[i],
                       col = piecol)
    }
  }
  if (!is.null(text))
    text(XX, YY, text, adj = adj, col = col, ...)
  if (!is.null(pch))
    points(
      XX + adj[1] - 0.5,
      YY + adj[2] - 0.5,
      pch = pch,
      col = col,
      bg = bg,
      ...
    )
}

# Additional functions
capwords = function(s, strict = TRUE) {
  cap <- function(s)
    paste(toupper(substring(s, 1, 1)),
          {
            s = substring(s, 2)
            if (strict)
              tolower(s)
            else
              s
          },
          sep = "", collapse = " ")
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
sanitize = function(str) {
  gsub('([#$%&~_\\^\\\\{}\\s\\(\\)])', '\\\\\\1', str, perl = TRUE)
}
getDistString = function (X, kwd) {
  topLeft = which(X == kwd, arr.ind = TRUE)
  if (length(topLeft) == 0) {
    dist = 'Delta' # default
    pars = '0'
  } else {
    dist = X[topLeft[1], topLeft[2] + 1]
    pars = X[topLeft[1], (topLeft[2] + 2):ncol(X)]
    pars = pars[!is.na(pars)]
  }
  string = paste0(dist, '(', paste0(pars, collapse = ','), ')')
}
getParams = function (X, kwd) {
  topLeft = which(X == kwd, arr.ind = TRUE)
  if (length(topLeft) == 0) {
    pars = NA
  } else {
    pars = X[topLeft[1], (topLeft[2] + 1):ncol(X)]
    if (any(!is.na(pars))) {
      pars = pars[!is.na(pars)]
    } else {
      pars = NA
    }
  }
}
sampleDistString = function (distString, sampleSize) {
  spl = unlist(strsplit(distString, '(', fixed = TRUE))
  if (spl[1] == 'Delta') {
    par1 = strsplit(spl[2], ',', fixed = TRUE)[1]
    par1 = sub(')', '', par1)
    sample = matrix(as.numeric(par1), ncol = 1, nrow = sampleSize)
  } else {
    sample = nds(sampleSize, distString)
  }
  return(sample)
}
nds = function(ns, dist) {
  command = paste("echo ", ns, " '", dist, "' | ./Rnested.x")
  # quotes around dist avoid shell interpretation
  tc = textConnection(system(command, intern = T))
  liste = scan(tc, quiet = TRUE)
  close(tc)
  nleaves = liste[1]
  nlast = nleaves * ns + 1
  nds = matrix(liste[2:nlast], ncol = nleaves, byrow = T)
  return(nds)
}
oneDist = function(ic, d, tags, tagged = FALSE) {
  dSub = d[[ic]]
  string = paste0(dSub$dist, '(')
  sSub = c()
  for (ip in 1:length(dSub$elem)) {
    sSub[ip] = paste0(
      dSub$mu[ip],
      ifelse(
        dSub$link[ip] == 0,
        ifelse(tagged,
               # paste0(":'", tags[dSub$elem[ip]], "'"),
               paste0("<", tags[dSub$elem[ip]], ">"),
               ''),
        paste0('*LINK/', dSub$link[ip], '/')
      ),
      collapse = '')
  }
  sig = dSub$sig[which(!is.na(dSub$sig))]
  return(paste0(
    string,
    paste0(sSub, collapse = ','),
    ifelse(length(sig) == 0,
           '',
           paste0(';', paste0(sig, collapse = ','))),
    ')'
  ))
}
oneNewick = function(ic, d, tags, tagged = TRUE) {
  dSub = d[[ic]]
  sSub = c()
  for (ip in 1:length(dSub$elem)) {
    sSub[ip] = paste0(
      ifelse(dSub$link[ip] == 0,
             tags[dSub$elem[ip]],
             paste0('LINK/', dSub$link[ip], '/')),
      collapse = '')
  }
  return(paste0("(", paste(sSub, collapse = ","), ")"))
}
oneEdgeTag = function(ic, d, tags) {
  dSub = d[[ic]]
  #   sig=dSub$sig[which(!is.na(dSub$sig))]
  mu = as.numeric(dSub$mu)
  sig = as.numeric(dSub$sig)
  sSub = c()
  for (ip in 1:length(dSub$elem)) {
    label = switch(
      dSub$dist,
      Dirg = {
        paste0(signif(mu[ip], 2), ' +/- ',
               signif(sig[ip], 2))
      },
      Mlgn = {
        paste0(signif(mu[ip], 2), ' */: ',
               signif(sig[ip], 2))
      },
      Diut = {
        paste0('[', signif(mu[ip], 2), '; ',
               signif(sig[ip], 2), ']')
      },
      Diri = {
        paste0(signif(mu[ip] / sum(mu), 2))
      },
      Dior = {
        paste0(signif((sum(mu) - mu[ip]) / sum(mu), 2))
      }
    )
    if (dSub$link[ip] == 0) {
      sSub[ip] = label
    } else {
      sSub[ip] = paste0(label, ',LINK/', dSub$link[ip], '/')
    }
  }
  return(paste0(sSub, collapse = ","))
}
oneNodeTag = function(ic, d, tags) {
  dSub = d[[ic]]
  #   sig=dSub$sig[which(!is.na(dSub$sig))]
  sSub = c()
  for (ip in 1:length(dSub$elem)) {
    if (dSub$link[ip] == 0) {
      sSub[ip] = dSub$dist[ip]
    } else {
      sSub[ip] = paste0(dSub$dist[ip], ',LINK/', dSub$link[ip], '/')
    }
  }
  return(paste0(sSub, collapse = ","))
}
getSpecies = function (chain) {
  species = unlist(strsplit(chain, ' + ', fixed = TRUE))
  species = as.vector(sapply(species, str_trim))
  return(species)
}
writeSample = function(iMC, dir, reac, tags, drawPars, drawBR, type) {
  # Generate csv file of a database draw
  dbOut = data.frame()
  reactants = getSpecies(reac)
  signBeta = 1
  if ('E' %in% reactants)
    signBeta = -1
  for (ip in 1:length(tags)) {
    prods = getSpecies(tags[ip])
    db1 = data.frame(
      R1 = reactants[1],
      R2 = reactants[2],
      R3 = reactants[3],
      P1 = prods[1],
      P2 = prods[2],
      P3 = prods[3],
      P4 = prods[4],
      a = drawPars['ALPHA'] * drawBR[ip],
      b = drawPars['BETA'] * signBeta,
      c = drawPars['GAMMA'],
      f = 1.0,
      g = 0.0,
      type = type
    )
    dbOut = rbind(dbOut, db1)
  }
  dbOutm <- within(dbOut, {
    b <- sprintf("%6.3f", b)
    a <- sprintf("%6.3e", a)
  })
  write.table(
    dbOutm,
    file = file.path(dir, paste0('run_', sprintf('%04i', iMC), '.csv')),
    quote = TRUE,
    sep = ';',
    na = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}
printBib = function(keys, bibList) {
  if (length(keys) != 0) {
    cat('<H2>References</H2><DL>\n')
    for (key in keys) {
      cat(paste0('<DT>[', key, ']</DT><DD>'))
      print(bibList[key], style = "html")
      cat('</DD>')
    }
    cat('</DL>')
  }
}
printBibKeys = function(keys) {
  if (any(!is.na(keys))) {
    refs = paste0(sort(keys), collapse = ',')
    paste0('[', refs, ']')
  } else {
    NA
  }
}
printRQ = function(comments) {
  if (!is.na(comments)) {
    cat('<H2>Comments</H2>\n')
    cat(paste0('<font color="red">', comments, '</font>\n'))
  }
}

# New BR representation in DB ####
getNewickFromTaggedDist = function(str) {
  # Get Newick string from tagged distribution
  newick = ''
  match = regexpr('<',str)
  while(match != -1) {
    sbstr = substr(str,1,1+match)
    newick = paste0(newick,gsub('[[:alnum:]*.;+<[:blank:]]{1}?','',sbstr))
    str  = substr(str,1+match,nchar(str))
    match = regexpr('>',str)
    sbstr = substr(str,1,match-1)
    newick = paste0(newick,sbstr)
    str  = substr(str,match+1,nchar(str))
    match = regexpr('<',str)
  }
  newick = paste0(newick,gsub('[[:alnum:]*.;+>,[:blank:]]{1}?','',str))
  newick = paste0(newick,';')
  return(gsub('-','',newick))
}
getDistFromTaggedDist = function(str1) {
  # Get dist string from tagged distribution
  gsub('<.*>{1}?','',str1)
}
getTagsFromTaggedDist = function(str1) {
  # Get Newick string from tagged distribution
  newick = getNewickFromTaggedDist(str1)
  unlist(strsplit(gsub(';','',gsub('\\(','',gsub('\\)','',newick))),','))
}
