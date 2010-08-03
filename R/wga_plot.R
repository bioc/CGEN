# History Mar 11 2008  Initial coding
#         Mar 13 2008  Add option to transform the locations 
#         Mar 18 2008  Check for errors in the beginning
#         Mar 27 2008  Add option to split the screen in 2.
#         Apr 08 2008  Add option for the maximum upper limit and
#                      a dashed line.
#         May 13 2008  Add code for UNIX.
#         Jul 25 2008  Add QQ.plot
#         Aug 21 2008  In QQ.plot, plot on -log10 scale
#         Sep 03 2008  Add plotting symbols and colors 
#         Oct 15 2008  Add options for QQ.plot
#         Nov 17 2008  Generalize chromosome.plot
#         Nov 20 2008  Remove missing values in chromosome, for the
#                      chromosome plot function.
#         Nov 25 2008  Add title and legend options in chromosome plot
#         Nov 25 2008  Fix bug with matching snp names in chromosome.plot
#         Mar 17 2009  Add forest.plot function
#                      Change to setDevice
#         Apr 11 2009  Add alternating colors to chromosome plot
#         Apr 27 2009  Allow alternating colors with multiple p-values
#                      in chromsome plot.
#         Apr 29 2009  Allow for character variables in qq, chromosome plot
#         Jul 15 2009  Add optional vec to QQ plot
#         Jul 23 2009  Remove library(rmeta) call in forestPlot
#         Jul 28 2009  Let alt.colors = 1 be default value in chromosome.plot
#                      Move forest.plot to wga_plot2.R
#         Oct 13 2009  Include add option in chromosome.plot
#                      Dashed line option
#                      Add options for padj and las in chromosome.plot
#         Oct 14 2009  Add getColors function
#         Nov 19 2009  Fix bug in chrm.plot.ylim
#         Mar 09 2010  Add functions set.plot and save.plot
#         Mar 15 2010  Add new option to chromosome plot for removing white space
#                      on edges of plot.

# Customize the plot.scan.gwaa function in GenABEL
myplot.scan.gwaa <- function (x, df=1, y.lower.lim=NULL, 
                    y.upper.lim=NULL, outfile=NULL,
                    type="jpeg", main=NULL, transform.loc=1,
                    P1df.pch=21, P2df.pch=24, Pc1df.pch=22,
                    xlab=NULL, legend=NULL, P1df.col="black",
                    P2df.col="red", Pc1df.col="green",
                    alternate.col=NULL, sub=NULL, hline=NULL) {
  
  # x              Input list with names "map", "P1df", 
  #                and possibly "P2df", "Pc1df", "chromosome"
  #                map is the vector of x-coordinates
  #                P1df, P2df, Pc1df are the vectors of y-coordinates
  # df             "1", "2", or "all" 
  #                The default is 1
  # y.lower.lim    Upper limit of p-values to plot or NULL
  # y.upper.lim    Lower limit of p-values or NULL
  # outfile        Output file
  # type           One of the following for the type of plot: 
  #                    ("wmf", "emf", "png", "jpeg", "jpg", "bmp",
  #                        "ps", "eps", "pdf")
  #                The default is "jpeg"
  # main           NULL or main plot title
  # transform.loc  0 or 1 to transform the snp locations.
  #                It is only done if the number of chromosomes is
  #                larger than 1.
  #                The default value is 1.
  # Pxxx.pch       19-25 plotting characters
  # xlab           x-axis label
  # legend         List with the names:
  #                where:   "bottomright", "bottom", "bottomleft", 
  #                         "left", "topleft", "top", "topright", 
  #                         "right" and "center". 
  #                         or specify a coordinate pair
  #                legend:   Character expression
  #                bty:      "o" or "n" type of box around the legend
  #                title:    title of legend 
  # Pxxx.col       Color used for the points
  # alternate.col  Vector of colors used to seperate chromosomes
  # splitScreen    0 or 1 to split the screen into two parts. Only valid
  #                for more than 1 chrmosome.
  # sub            Subtitle (at the bottom of plot)
  #                The default is NULL
  # hline          y value to draw a dashed horizontal line
  #                The default is NULL

    df <- as.character(df)
    if (!(any(c("1", "2", "Pc1df", "all") == df))) 
        stop("df parameter must be 1, 2, \"Pc1df\" or \"all\"")

    # Check vectors
    if (is.null(x$map)) stop("ERROR: map must be specified")
    n <- length(x$map)
    if (is.null(x$P1df)) {
      x$chromosome <- factor(rep(1, n))
    } else if (n != length(x$P1df)) {
      stop("length of map and P1df not equal!")
    }
    if (is.null(x$chromosome)) {
      x$chromosome <- factor(rep(1, n))
    } else if (n != length(x$chromosome)) {
      stop("length of map and chromosome not equal!")
    } else {
      if (!is.factor(x$chromosome)) x$chromosome <- factor(x$chromosome)
    }
    if (is.null(x$P2df)) {
      x$P2df <- factor(rep(NA, n))
    } else if (n != length(x$P2df)) {
      stop("length of map and P2df not equal!")
    }
    if (is.null(x$Pc1df)) {
      x$Pc1df <- factor(rep(NA, n))
    } else if (n != length(x$Pc1df)) {
      stop("length of map and Pc1df not equal!")
    }

    if (any(names(x) == "Pgw1df")) { 
      if (!is.null(x$Pgw1df)) {
        x$P1df <- x$Pgw1df
        x$P2df <- x$Pgw2df
      }
    }

    # Save vectors
    formula    <- x$formula
    map        <- x$map
    chromosome <- x$chromosome
    P1df       <- x$P1df
    P2df       <- x$P2df
    Pc1df      <- x$Pc1df

    # Free memory
    rm(x)
    temp <- gc()

    if (df == "1") {
      Pv  <- P1df
      pch <- P1df.pch
      col <- P1df.col
      rm(P2df, Pc1df)
      temp <- gc()
    }
    else if (df == "2") {
      Pv  <- P2df
      pch <- P2df.pch
      col <- P2df.col
      rm(P1df, Pc1df)
      temp <- gc()
    }
    else if (df == "Pc1df") {
      Pv  <- Pc1df
      pch <- Pc1df.pch
      col <- Pc1df.col
      rm(P1df, P2df)
      temp <- gc()
    } else {
      Pv  <- P1df
      pch <- P1df.pch
      col <- P1df.col
    }
    Pv       <- replace(Pv, (Pv <= 0), 1e-16)
    mLog10Pv <- -log10(Pv)

    if (df == "all") {
      P2df  <- replace(P2df, (P2df <= 0), 1e-16)
      Pc1df <- replace(Pc1df, (Pc1df <= 0), 1e-16)
      max1y <- max(-log10(Pv), na.rm = TRUE)
      max2y <- max(-log10(P2df), na.rm = TRUE)
      max3y <- max(-log10(Pc1df), na.rm = TRUE)
      maxy  <- max(c(max1y, max2y, max3y), na.rm = TRUE)
    }
    else {
      maxy  <- max(mLog10Pv, na.rm=TRUE)
    }

    # Get ylim if not passed in
    if (is.null(y.lower.lim)) {
      y.lower.lim <- 0
    } else {
      y.lower.lim <- -log10(y.lower.lim)
    }

    if (!is.null(y.upper.lim)) {
      y.upper.lim <- -log10(y.upper.lim) + 0.5 
    } else {
      y.upper.lim <- maxy
    }

    ylim <- c(y.lower.lim, y.upper.lim)

    # Get the main title
    if (is.null(main)) main <- formula

    chrms <- levels(chromosome)

    nchrm <- length(chrms)
    j     <- 1
    add0  <- max(as.integer(chrms), na.rm=TRUE)

    if (is.null(xlab)) {
      if (nchrm > 1) {
        xlab <- "Chromosome"
      } else {
        xlab <- "Map Position"
      }
    }

    # Remove Pv
    rm(Pv)
    temp <- gc()

    if (nchrm > 1) {
      # Transform locations
      if (transform.loc) {
        # Scale for each chromosome
        chpos  <- rep(NA, nchrm)
        addpos <- rep(NA, nchrm)
        for (i in 1:nchrm) {
          temp <- (chromosome == chrms[i])

          # Scale between 0 and 1 and translate
          max <- max(map[temp])
          add <- as.integer(chrms[i])
          if (is.na(add)) {
            add <- add0 + j
            j   <- j + 1
          }
          map[temp] <- map[temp]/max + add

          # Save these numbers add, we need them later
          addpos[i] <- add

          # Get the new mean
          chpos[i]  <- mean(map[temp])
        }
      } # END: if (transform.loc)

        if (df == "all") 
            plot(map, mLog10Pv, xlab=xlab,
                ylab=expression(-log[10](P-value)), axes=FALSE, 
                ylim=ylim, pch=pch, col=col)
        else plot(map, mLog10Pv, xlab=xlab,
            ylab=expression(-log[10](P-value)), axes=FALSE, ylim=ylim,
                pch=pch, col=col)
        #mxlog <- floor(max(mLog10Pv))
        mxlog <- floor(y.upper.lim)

        if (mxlog == 0) mxlog <- max(mLog10Pv)
        if (mxlog < 1) {
          axis(2, at=c(y.lower.lim, mxlog))
        } else {
          axis(2, at=floor(y.lower.lim:mxlog))
        }
        axis(1, at=chpos, labels=chrms)
        box()
    }
    else {
      if (df == "all") { 
        plot(map, mLog10Pv, xlab=xlab, 
             ylab = expression(-log[10](P-value)), ylim=ylim,
             col=col, pch=pch)
      } else {
        plot(map, mLog10Pv, xlab=xlab, 
             ylab = expression(-log[10](P-value)), ylim=ylim,
             col=col, pch=pch)
      }
    }

    # Add the title
    title(main=main, sub=sub)

    if (!is.null(alternate.col)) {
      alternate.col <- rep(alternate.col, 
                        times=ceiling(nchrm/length(alternate.col))) 
   
      if (df == "all") {
        pch  <- P1df.pch
        Pv2  <- -log10(P2df)
        pch2 <- P2df.pch
      } else if (df == "1") {
        pch  <- P1df.pch
      } else {
        pch  <- P2df.pch
      }

      # Sort the addpos vector
      temp <- sort(addpos, index.return=TRUE)$ix
      chrms <- chrms[temp]

      for (i in 1:nchrm) {
        temp <- (chromosome == chrms[i])
        points(map[temp], mLog10Pv[temp], pch=pch, col=alternate.col[i])

        if (df == "all") {
          points(map[temp], Pv2[temp], pch=pch2, col=alternate.col[i])
        }
      }
    } else {
      if (df == "all") {
        points(map, -log10(P2df), col=P2df.col, pch=P2df.pch)
        points(map, -log10(Pc1df), col=Pc1df.col, pch=Pc1df.pch)
      }
    }
    
    # Add a legend
    if ((!is.null(legend)) && (df == "all")) {
      pch <- c(P1df.pch, P2df.pch, Pc1df.pch)
      if (!is.null(alternate.col)) {
        col <- "black"
      } else {
        col <- c(P1df.col, P2df.col, Pc1df.col)
      }
      legend(x=legend$where, legend=legend$legend,
             bty=legend$bty, title=legend$title,
             pch=pch, col=col)
    }

    # Draw horizontal dashed (lty=2) line
    if (!is.null(hline)) {
      hline <- -log10(hline)
      abline(h=hline, lty=2)
    } 

    # Save the plot
    if (!is.null(outfile)) savePlot(filename=outfile, type=type)

    0


} # END: myplot.scan.gwaa

myplot <- function(x, df=1, y.lower.lim=NULL, y.upper.lim=NULL,
                    outfile=NULL,
                    type="jpeg", main=NULL, transform.loc=1,
                    P1df.pch=21, P2df.pch=24, Pc1df.pch=22,
                    xlab=NULL, legend=NULL, P1df.col="black",
                    P2df.col="red", Pc1df.col="green",
                    alternate.col=NULL, splitScreen=0, sub=NULL,
                    hline=NULL) {

  # splitScreen      0 or 1 to split the plot in 2 parts
  #                  The default is 0

  # Set the output device
  setDevice(outfile, list(type=type, which=0)) 

  x <- as.list(x)
  chrms <- levels(x$chromosome)
  nchrm <- length(chrms)
  if (nchrm <= 1) splitScreen <- 0

  if (splitScreen) {
    half <- floor(nchrm/2)
    c1 <- as.numeric(chrms)
    c1 <- sort(c1)
    chrm1 <- c1[1:half]

    xx   <- x
    temp <- x$chromosome %in% chrm1

    x$chromosome <- factor(x$chromosome[temp])

    x$map <- x$map[temp]    
    if (!is.null(x$P1df)) {
      x$P1df <- x$P1df[temp]
    }
    if (!is.null(x$P2df)) {
      x$P2df <- x$P2df[temp]
    }
    if (!is.null(x$Pc1df)) {
      x$Pc1df <- x$Pc1df[temp]
    }

    split.screen(c(2,1))
    screen(1)

    ret <- myplot.scan.gwaa(x, df=df, y.lower.lim=y.lower.lim, 
            y.upper.lim=y.upper.lim, outfile=NULL,
            type=type, main=main, transform.loc=transform.loc,
            P1df.pch=P1df.pch, P2df.pch=P2df.pch, Pc1df.pch=Pc1df.pch,
            xlab=xlab, legend=legend, P1df.col=P1df.col,
            P2df.col=P2df.col, Pc1df.col=Pc1df.col,
            alternate.col=alternate.col, sub=sub, hline=hline)

    temp <- as.logical(1-temp)
    x    <- xx

    x$chromosome <- factor(x$chromosome[temp])
    x$map <- x$map[temp]    
    if (!is.null(x$P1df)) {
      x$P1df <- x$P1df[temp]
    }
    if (!is.null(x$P2df)) {
      x$P2df <- x$P2df[temp]
    }
    if (!is.null(x$Pc1df)) {
      x$Pc1df <- x$Pc1df[temp]
    }

    screen(2)

    ret <- myplot.scan.gwaa(x, df=df, y.lower.lim=y.lower.lim, 
            y.upper.lim=y.upper.lim, outfile=NULL,
            type=type, main=main, transform.loc=transform.loc,
            P1df.pch=P1df.pch, P2df.pch=P2df.pch, Pc1df.pch=Pc1df.pch,
            xlab=xlab, legend=legend, P1df.col=P1df.col,
            P2df.col=P2df.col, Pc1df.col=Pc1df.col,
            alternate.col=alternate.col, sub=sub, hline=hline)


  
  } # END: if (splitScreen)
  else {
    ret <- myplot.scan.gwaa(x, df=df, y.lower.lim=y.lower.lim, 
            y.upper.lim=y.upper.lim, outfile=NULL,
            type=type, main=main, transform.loc=transform.loc,
            P1df.pch=P1df.pch, P2df.pch=P2df.pch, Pc1df.pch=Pc1df.pch,
            xlab=xlab, legend=legend, P1df.col=P1df.col,
            P2df.col=P2df.col, Pc1df.col=Pc1df.col,
            alternate.col=alternate.col, sub=sub, hline=hline)
  }


  # Save the plot
  setDevice(outfile, list(type=type, which=1)) 

  0

} # END: myplot

# Function to set a graphics device
setDevice <- function(file, op=NULL) {

  # file     NULL or file to save plot
  ###############################################################
  # op
  #  type    "jpeg", "ps" or "pdf"
  #          The default is "jpeg"
  #  which   0 or 1  0=before plot, 1 = after plot
  #          The default is 0

  if (is.null(file)) return(NULL)

  op   <- default.list(op, c("type", "which"), list("jpeg", 0))
  type <- tolower(op$type)

  winFlag <- (.Platform$OS.type == "windows")

  if (op$which == 0) {
    if (winFlag) {
      return(0)
    }
    if (type == "ps") {
      postscript(file)
    } else if (type == "pdf") {
      pdf(file)
    } else {
      jpeg(file)
    } 
  } else {
    if (!winFlag) {
      graphics.off()
    } else {
      savePlot(filename=file, type=type)
    }
  }

  0

} # END: setDevice

# Function to create a QQ plot on a -log10 scale.
QQ.plot <- function(pvals, op=NULL) {

  # pvals      Vector of p-values
  ###############################################
  # op         List with the optional fields:
  #   title    Title for the plot
  #            The default is "QQ Plot"
  #   ylim     NULL or vector of length 2.
  #            ylim=c(0, 10) will cause the y-axis range to be between
  #            10^{-0} and 10^{-10}
  #            The default is NULL
  #   outfile  File to save the plot or NULL.
  #            The default is NULL
  #   type     "jpeg", "ps", or "pdf"
  #            The default is "jpeg"
  ###############################################

  op <- default.list(op, c("title", "type", "color"), 
                     list("QQ Plot", "jpeg", "blue"))
  op$type <- tolower(op$type)

  pvals <- as.numeric(pvals)

  # Remove non-finite values
  pvals <- pvals[is.finite(pvals)]

  pvals <- sort(pvals)
  n     <- length(pvals)
  ex    <- (1:n)/(n+1)

  # Plot on -log10 scale
  pvals[pvals < 1e-16] <- 1e-16
  pvals <- -log10(pvals)
  ex    <- -log10(ex)

  xlabel <- expression(-log[10]("Expected P-values"))
  ylabel <- expression(-log[10]("Observed P-values"))

  temp <- setDevice(op$outfile, op=list(type=op$type, which=0)) 

  plot(ex, pvals, axes=FALSE, type="p", xlab=xlabel, ylab=ylabel,
        lwd=1, col=op$color[1], ylim=op$ylim)

  maxy <- ceiling(max(pvals))
  maxx <- ceiling(max(ex))
  m    <- max(maxy, maxx)
  my   <- m
  if (!is.null(op$ylim)) my <- max(m, op$ylim)
  ypos <- 0:my
  xpos <- 0:m

  axis(1, at=xpos)
  axis(2, at=ypos)
  box()
  title(main=op$title, col.main="blue", cex.main=1)

  # Add a diagonal line
  abline(a=0, b=1, col="black", lwd=2)

  # See if other p-values need to be plotted
  pvals <- getListName(op, "pvalues")
  if (!is.null(pvals)) {
    pvals <- as.numeric(pvals)
    pvals <- sort(pvals)
    n     <- length(pvals)
    ex    <- (1:n)/(n+1)
    pvals[pvals < 1e-16] <- 1e-16
    pvals <- -log10(pvals)
    ex    <- -log10(ex)
    points(ex, pvals, type="p", lwd=1, col=op$color[2])
  }

  temp <- setDevice(op$outfile, op=list(type=op$type, which=1)) 

  0 

} # END: QQ.plot

# Function to transform locations for a chromosome plot
transform.loc <- function(chromosome, chrms, map, addToEnd=0) {

  nchrm <- length(chrms)

  # Get the maximum value and length of each chrm
  clen <- rep(NA, times=nchrm)
  cmax <- clen
  for (i in 1:nchrm) {
    temp <- map[(chromosome == chrms[i])]

    # Get the max and min values
    cmax[i] <- max(temp, na.rm=TRUE)
    clen[i] <- cmax[i] - min(temp, na.rm=TRUE) 
  }

  # Remove problem values
  xx <- !is.finite(cmax)
  if (any(xx)) {
    for (i in 1:nchrm) {
      if (xx[i]) {
        temp <- !(chromosome == chrms[i])
        chromosome <- chromosome[temp]
        map        <- map[temp]
      }
    }
    chrms <- chrms[!xx]
  }

  # Scale the lengths
  clen  <- clen/max(clen)
  scale <- clen/cmax 

  # Scale for each chromosome
  chpos <- rep(NA, nchrm)
  add   <- 0
  cmin  <- rep(NA, nchrm)

  for (i in 1:nchrm) {
    temp <- (chromosome == chrms[i])

    # Scale between 0 and 1 and translate
    map[temp] <- scale[i]*map[temp] + add

    # Get min and max
    cmin[i] <- min(map[temp], na.rm=TRUE)
    cmax[i] <- max(map[temp], na.rm=TRUE)

    # Get the new mean
    chpos[i] <- (cmin[i] + cmax[i])/2

    # Update add
    add <- cmax[i] + addToEnd
  }

  list(map=map, chpos=chpos, cmin=cmin, cmax=cmax)

} # END: transform.loc

# Function to create a chromosome plot
chromosome.plot <- function(infile, plot.vars, locusMap.list, op=NULL) {

  # infile           Input data set containing the p-values.
  #                  infile can be the path to the file or a data
  #                  frame containing the pvalues and the locus map
  #                  variables.
  #                  No default.
  # plot.vars        Variables in infile to plot
  # locusMap.list    List containing info about the file containing, 
  #                  snp, chromosome, and location.
  ###################################################################
  # op               List of options:
  #  splitScreen     0 or 1 to split the plot into 2 parts
  #                  The default is 1.
  #  snp.var         The snp variable name in infile
  #                  The default is "SNP".
  #  yaxis.range     Range for the y-axis. Should be on the original scale.
  #                  Ex: c(1e-12, 1e-3)
  #                  The default is NULL.
  #  title           NULL or title of plot
  #                  The default is NULL
  #  legend          0 or 1 for a legend
  #                  The default is 1
  #  legend.names    The default is plot.vars
  #  legend.horiz    TRUE or FALSE for a horizontal legend
  #                  The default is 1
  #  legend.where    "bottomright", "bottom", "bottomleft", 
  #                         "left", "topleft", "top", "topright", 
  #                         "right" and "center". 
  #                  The default is "top"
  #  cex.axis        X-axis label size
  #                  The default is 0.75
  #  colors          Vector of colors to use
  #                  The default is NULL
  #  alt.colors      0 or 1 to alternate colors in plot
  #                  The default is 1
  #  pch             Vector of plotting symbols
  #                  The default is 21 ("circles")
  #                  19 = solid circle, 20 = bullet
  #  add             A number to add spacing between the chromsomes
  #                  The default is 0.
  #  x.padj          X-axis padj option
  #                  The default depends on x.las
  #  x.las           0-3 for x-axis labels 
  #                  0=parallel, 1=horizontal, 2=perpendicular, 3=vertical
  #                  The default is 0
  #  xlim.add        Vector of length 1 or 2 for adding(subtracting) a 
  #                  value to xlim
  #                  The default is c(0, 0)
  #####################################################################
  #  hline           NULL or list specifying a horizontal line
  #    h             y value
  #    lty
  #                  The default is NULL
  #####################################################################

  op      <- default.list(op, 
            c("splitScreen", "snp.var", "legend", "cex.axis", "alt.colors",
              "legend.horiz", "legend.where", "add", "x.las", "xlim.add"), 
             list(1, "SNP", 1, 0.75, 1, "TRUE", "top", 0, 0, c(0, 0)))
  subset  <- getListName(op, "subset")
  subFlag <- !is.null(subset)

  plot.vars <- unique(plot.vars)
  nvars     <- length(plot.vars)

  dfFlag  <- is.data.frame(infile) | is.matrix(infile)

  if (!dfFlag) {
    # Check the locusMap list
    locusMap.list <- check.locusMap.list(locusMap.list)

    # Read in the locus map data
    snp  <- NULL
    chrm <- NULL
    map  <- NULL
    for (file in locusMap.list$file) {
      temp <- paste(locusMap.list$dir, file, sep="")
      data <- getLocusMap(temp, locusMap.list)
      snp  <- c(snp, data$snp)
      chrm <- c(chrm, data$chrm)
      map  <- c(map, data$loc)
    }
    rm(data)
    temp <- gc(verbose=FALSE)

  } else {
    infile <- unfactor.all(infile)
    snp    <- infile[, locusMap.list$snp.var]
    chrm   <- as.character(infile[, locusMap.list$chrm.var])
    map    <- as.numeric(infile[, locusMap.list$loc.var])
  }

  if (subFlag) {
    temp <- chrm %in% subset
    chrm <- chrm[temp]
    snp  <- snp[temp]
    map  <- map[temp]
  }

  if (!dfFlag) {
    # Read in the p-values
    tlist <- list(file=infile, file.type=3, header=1, delimiter="\t")
    data  <- getColumns(tlist, c(op$snp.var, plot.vars), temp.list=NULL)

    # Match the snp names
    snp2 <- data[[op$snp.var]]
    data[[op$snp.var]] <- NULL
    rm(tlist)
  } else {
    snp2 <- unfactor(infile[, op$snp.var])
    data <- list()
    for (var in plot.vars) data[[var]] <- as.numeric(infile[, var])
    rm(infile)
    temp <- gc(verbose=FALSE)
  }

  # Get the correct subset for data
  temp  <- snp2 %in% snp
  snp2  <- snp2[temp]
  for (var in plot.vars) {
    data[[var]] <- as.numeric(data[[var]][temp])
  }

  # Match chrm and location to the data
  temp  <- match(snp2, snp)
  if (any(is.na(temp))) stop("ERROR: matching SNP names")
  chrm  <- chrm[temp]
  map   <- map[temp]

  maxp <- rep(NA, times=nvars)
  # Transform p-values to a -log10 scale
  for (i in 1:nvars) { 
    temp <- data[[i]] < 1e-30
    data[[i]][temp] <- 1e-30
    data[[i]]       <- -log10(data[[i]])
    maxp[i]         <- max(data[[i]], na.rm=TRUE)
  } 

  # Remove values
  temp <- (is.finite(map)) & (!is.na(chrm))
  for (i in 1:nvars) temp <- (temp & is.finite(data[[i]]))
  map  <- map[temp]
  chrm <- chrm[temp]
  for (i in 1:nvars) data[[i]] <- data[[i]][temp]

  # y-axis limits
  temp <- getListName(op, "yaxis.range")
  ylim <- chrm.plot.ylim(maxp, ylim=temp)

  rm(maxp, var, snp, snp2)
  temp <- gc(verbose=FALSE)
  uchrm <- getUniqueChrm(chrm)
  nchrm <- length(uchrm)

  if (nchrm <= 1) {
    op$splitScreen <- 0
    transLoc       <- 0
  } else {
    transLoc       <- 1
  }

  # axis labels
  if (nchrm > 1) {
    xlab <- "Chromosome"
  } else {
    xlab <- "Map Position"
  }
  ylab <- expression(-log[10](P-value))

  # Legend
  nclrs <- nvars
  alt.colors <- op$alt.colors
  if (alt.colors) nclrs <- nchrm
  colors <- getListName(op, "colors")
  colors <- chrm.plot.colors(nclrs, colors=colors)
  pch    <- chrm.plot.pch(nvars, pch=op$pch)
  if (op$legend) {
    names <- plot.vars
    temp <- getListName(op, "legend.names")
    if (!is.null(temp)) names <- temp
    legend <- chrm.plot.legend(names, colors, pch, alt.colors,
               horiz=op$legend.horiz, where=op$legend.where)
  } else {
    legend <- NULL
  }

  hline  <- getListName(op, "hline")
  x.padj <- getListName(op, "x.padj") 

  if (op$splitScreen) {

    half  <- floor(nchrm/2)
    chrm1 <- uchrm[1:half]

    # Get the ids
    temp  <- chrm %in% chrm1
    
    # Top plot
    split.screen(c(2,1))
    screen(1)

    ret <- chrm.plot.main(data, map, chrm, ids=temp, ylim=ylim[1:2],
           xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
           colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
           alt.colors=alt.colors, add=op$add, hline=hline,
           x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)

    # Get the ids for the bottom plot
    temp <- as.logical(1-temp)

    screen(2)
    ret <- chrm.plot.main(data, map, chrm, ids=temp, ylim=ylim[3:4],
            xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors=alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)

  } # END: if (splitScreen)
  else {
    ret <- chrm.plot.main(data, map, chrm, ids=NULL, ylim=ylim[1:2],
              xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)
  }

} # END: chromosome.plot

# Function to produce a plot
chrm.plot.main <- function(pvals, map, chrm, ids=NULL, ylim=c(0, 8),
                  xlab=NULL, ylab=NULL, transLoc=1, legend=NULL,
                  colors=NULL, pch=NULL, title=NULL, cex.axis=0.75,
                  alt.colors=0, add=0, hline=NULL, x.las=0, x.padj=-1.0,
                  xlim.add=c(0,0)) {

  # pvals   List
 
  if (is.null(ids)) ids <- 1:length(map)
  if (is.null(x.padj)) {
    if (x.las %in% c(0, 1)) {
      x.padj <- -1.0
    } else {
      x.padj <- 0.5
    }
  }

  # Subset the data
  chrm  <- chrm[ids]
  map   <- map[ids]
  vars  <- names(pvals)
  nvars <- length(vars)

  for (i in 1:nvars) pvals[[i]] <- pvals[[i]][ids]

  rm(ids)
  temp <- gc(verbose=FALSE)
  
  uchrms <- getUniqueChrm(chrm)
  nchrms <- length(uchrms)

  if (transLoc) {
    temp   <- transform.loc(chrm, uchrms, map, addToEnd=add)
    map    <- temp$map
    chpos  <- temp$chpos
    cmin   <- temp$cmin
    cmax   <- temp$cmax
  }

  nclrs <- nvars
  if (alt.colors) nclrs <- nchrms

  # Get the colors and plotting characters
  col  <- chrm.plot.colors(nclrs, colors=colors)
  pch  <- chrm.plot.pch(nvars, pch=pch)
  if (length(xlim.add) == 1) xlim.add <- rep(xlim.add, times=2)
  xlim <- c(min(map, na.rm=TRUE)+xlim.add[1], max(map, na.rm=TRUE)+xlim.add[2])

  # Initialize the plot
  if (transLoc) {
    
    plot(map, pvals[[1]], xlab=xlab,ylab=ylab, axes=FALSE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
    mxlog <- floor(ylim[2])
    #if (mxlog == 0) mxlog <- maxp[1]
    if (mxlog < 1) {
      axis(2, at=c(ylim[1], mxlog))
    } else {
      axis(2, at=floor(ylim[1]:mxlog))
    }

    # tcl is for tick mark length
    # padj is for moving the labels close/farther from the axis
    axis(1, at=chpos, labels=uchrms, tcl=-0.5, padj=x.padj, las=x.las,
         cex.axis=cex.axis, font=2)
    box()   

  } else {
    plot(map, pvals[[1]], xlab=xlab,ylab=ylab, axes=TRUE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
  }

  # Plot the other vars
  if (nvars > 1) {
    for (i in 2:nvars) {
      points(map, pvals[[i]], col=col[i], pch=pch[i])
    }
  }

  # Alternate colors
  if (alt.colors) {
    for (j in 1:nvars) {
      for (i in 1:nchrms) {
        temp <- (chrm == uchrms[i])
        points(map[temp], pvals[[j]][temp], pch=pch[j], col=col[i])
      }
    }
  }

  # Add line segments (tick marks)
  if ((transLoc) & (!alt.colors)) addLineSegments(cmin, cmax, ylim) 

  # Add a legend
  if (!is.null(legend)) {
    legend(x=legend$x, y=legend$y, legend=legend$legend,
           bty=legend$bty, title=legend$title,
           pch=legend$pch, col=legend$col, cex=legend$cex, 
           horiz=legend$horiz)
  }

  # Add title 
  if (!is.null(title)) title(title)

  if (!is.null(hline)) {
    hline <- default.list(hline, c("h", "lty"), list(7, 2))
    abline(h=hline$h, lty=hline$lty)
  }

} # END: chrm.plot.main

# Function to get the unique chromosomes
getUniqueChrm <- function(chrm) {

  u    <- unique(chrm)
  n    <- as.numeric(u)
  nas  <- is.na(n)
  ch   <- u[nas]
  s    <- sort(n)
  c(s, ch)

} # END: getUniqueChrm

# Function to return a vector of colors
chrm.plot.colors <- function(n, colors=NULL) {

  if (is.null(colors)) {
    col    <- c("blue", "pink", "green", "red", "black", "orange",
                "yellow", "gray", "brown", "turquoise", "gold",
                "violet", "olivedrab", "skyblue", "purple",
                "maroon")
    col <- c(col, col)
    #col <- rainbow(n)
  } else {
    col <- rep(colors, times=n)
  }
  col    <- col[1:n]
  col

} # END: chrm.plot.colors

# Function to return a vector of plotting characters
chrm.plot.pch <- function(n, pch=NULL) {

  if (is.null(pch)) {
    pch <- rep(21, times=n)
  } else {
    
    pch <- rep(pch, times=n)
  }
  pch <- pch[1:n]
  pch

} # END: chrm.plot.pch

# Function to return the legend list
chrm.plot.legend <- function(plot.vars, colors, pch, alt.colors,
                             horiz=TRUE, where="top") {

  nvars  <- length(plot.vars)
  col    <- chrm.plot.colors(nvars, colors=colors)
  if ((alt.colors) && (nvars > 1)) col <- rep("black", times=nvars)
  title  <- NULL
  x      <- where
  y      <- NULL
  legend <- plot.vars
  cex    <- 0.7
  horiz  <- horiz

  list(x=x, y=y, bty="n", pch=pch, col=col, title=title, 
       legend=legend, cex=cex, horiz=horiz)  

} # END: chrm.plot.legend

# Function to add line segments to a graph
addLineSegments <- function(cmin, cmax, ylim) {

  n    <- length(cmin) + 1
  temp <- cmin 
  temp[1] <- cmin[1] - 0.035
  temp[n] <- cmax[n-1] + 0.035
  if (n > 2) {
    for (i in 2:(n-1)) {
      temp[i] <- (cmax[i-1] + cmin[i])/2
    }
  }
  temp0 <- rep.int(-1, times=n)
  temp1 <- rep.int(ylim[1]+0.1, times=n)
  segments(temp, temp0, temp, temp1, lwd=1, font=2)

} # END: addLineSegments

# Function to return the y-axis limits
chrm.plot.ylim <- function(maxp, ylim=NULL) {

  # maxp is on a -log10 scale
  temp <- max(maxp) + 1
  y0   <- c(0, temp, 0, temp)
  if (is.null(ylim)) return(y0)
  n <- length(ylim)
  if (!(n %in% c(2, 4))) return(y0)
  
  for (i in 1:n) {
    if ((ylim[i] > 0) && (ylim[i] <= 1)) ylim[i] <- -log10(ylim[i])
  }

  # check ylim
  if (ylim[1] > ylim[2]) {
    temp    <- ylim[1]
    ylim[1] <- ylim[2]
    ylim[2] <- temp
  }
  if (n == 4) {
    if (ylim[3] > ylim[4]) {
      temp    <- ylim[3]
      ylim[3] <- ylim[4]
      ylim[4] <- temp
    }
  }

  if (n == 2) ylim <- c(ylim, ylim)

  ylim

} # END: chrm.plot.ylim

# Function to return colors
getColors <- function(colVec, op=NULL) {

 # op      List with names
 #  ncolors  (Maximum) number of colors to return
 #           Default is NULL
 #  plot     0 or 1
 #  print    0 or 1
 #  seed     Default is -1
 #  exclude  Character vector of colors to exclude

 op <- default.list(op, c("plot", "print", "seed", "sample"), list(0, 0, -1, 0))

 all <- colors()
 cls <- NULL
 for (cc in colVec) {
   temp <- grep(cc, all, fixed=TRUE)
   if (length(temp)) cls <- c(cls, all[temp])
 }
 
 rm(all)
 gc()
 exclude <- getListName(op, "exclude")
 if (!is.null(exclude)) {
   for (cc in exclude) {
     temp <- grep(cc, cls, fixed=TRUE)
     if (length(temp)) cls <- cls[-temp]
   }
 } 

 seed <- op$seed
 if (seed > 0) set.seed(seed)
 n <- length(cls)
 ncolors <- getListName(op, "ncolors")
 if (!is.null(ncolors)) {
   # If n != ncolors, sample from the colors
   if (n != ncolors) {
     if (n > ncolors) {
       cls <- sample(cls, ncolors, replace=FALSE)
     } else {
       temp <- ncolors - n
       if (temp > n) {
         temp <- sample(cls, temp, replace=TRUE)
       } else {
         temp <- sample(cls, temp, replace=FALSE)
       }
       cls <- c(cls, temp) 
     }
   }
 } else {
   ncolors <- n
 } 

 if (op$sample) cls <- sample(cls, ncolors, replace=FALSE)
 if (op$print) print(cls)
 if (op$plot) pie(rep.int(1, ncolors), col=cls)

 cls

} # END: getColors

# Function to call before plotting
set.plot <- function(op=NULL) {

  winFlag <- (.Platform$OS.type == "windows")
  if (winFlag) return(0)

  op <- default.list(op, c("out", "type", "res"),
           list("ERROR", "jpeg", 100), error=c(1, 0, 0))

  host <- callOS("hostname", intern=TRUE)
  temp <- op[["R_GSCMD", exact=TRUE]]
  if (is.null(temp)) {
    if (host == "biowulf.nih.gov") {
      temp <- "/data/wheelerb/nilanjan/wga/software/ghostscript-8.64/bin/gs" 
    } else {
      temp <- "/data/software/ghostscript-8.64/bin/gs"
    }
  }

  # Set the R_GSCMD environment variable for plots
  Sys.setenv(R_GSCMD=temp)

  bitmap(op$out, type=op$type, res=op$res)

  0 

} # END: set.plot

# Function to call after plotting
save.plot <- function(op=NULL) {

  winFlag <- (.Platform$OS.type == "windows")
  if (!winFlag) {
    #dev.off()
    return(0)
  }

  op <- default.list(op, c("out", "type"),
           list("ERROR", "jpeg"), error=c(1, 0))

  savePlot(filename=op$out, type=op$type)

  0 

} # END: save.plot
