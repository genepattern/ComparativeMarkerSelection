########################################################################
#                              MISCELLANEA                             #
########################################################################

VERBOSE <- function( v, ... )
{
  if ( v ) cat( ... )
}
my.nlevels <- function( x )
{
  return( length(my.levels(x)) )
}
my.levels <- function( x, sort=T )
{
  #if ( !is.numeric(x) ) stop( "x must be numeric" )
  return( if (sort) sort(unique(as.vector(x))) else unique(as.vector(x)) )
}
########################################################################
#                      READING RES/GCT/CLS FILES                       #
########################################################################
#
read.res <- function( file, nrow=T, sep="\t", verbose=F, complete=F, postfix=".dup" )
{
  calls <- desc <- scale <- NULL

  skip <- if (nrow) 3 else 2
  ncols <- (count.fields(file,sep=sep,comment.char="")[skip+1]-2)/2

  VERBOSE( verbose, "\tReading experiment names .. " )
  x.colnames <- read.table(file,sep=sep,fill=T,nrows=1,
                           check.names=F,comment.char="",quote="")[,-c(1,2)]
  VERBOSE( verbose, "done.\n" )

  VERBOSE( verbose, "\tReading signal .. " )
  x <- read.delim(file,header=F, sep=sep, fill=T, row.names=2, skip=skip,
                  check.names=F, comment.char="", quote="",
                  colClasses=c("character","character",rep(c("numeric","character"),ncols)))[,-1]
  VERBOSE( verbose, "done, [", nrow(x), " x ", ncol(x)/2, "] data entries.\n", sep="" )
  levs <- seq( 1,ncol(x),2 )

  if ( complete ) {
    VERBOSE( verbose, "\tReading calls .. " )
    calls <- as.matrix(x[,-levs])
    colnames(calls) <- sapply(x.colnames[levs],as.character)
    desc  <- as.character(read.delim(file,header=F,sep=sep,fill=T, skip=skip,check.names=F,
                                     comment.char="", quote="")[,1])
    scale <-
      sapply(read.table(file,sep=sep,fill=T,nrows=1,skip=1,
                        check.names=F,comment.char="",quote="")[1,levs+2],
             as.character)
    VERBOSE( verbose, "done.\n" )
  }
  x <- x[,levs]
  levs.na <- apply(x,2,function(z){ all(is.na(z)) })
  x <- x[,!(levs.na)]
  colnames(x) <- sapply(x.colnames[levs],as.character)
  x <- as.matrix(x)

  if ( complete ) {
    return (list(data=x,calls=calls,desc=desc,scale=scale))
  }
  else {
    return( x )
  }
}
read.gct <- function( file, complete=F, verbose=F )
{
  VERBOSE( verbose, "\tReading signal .. " )
  x <- read.delim(file,sep="\t",fill=T,skip=2, header=T,check.names=F,comment.char="")
  rownames(x) <- x[,1]
  desc <- as.character(x[,2])
  x <- as.matrix(x[,-c(1,2)])
  VERBOSE( verbose, "done, ", nrow(x), " genes x ", ncol(x), " experiments.\n", sep="" )

  if ( complete )
    return( list( data=x, desc=desc ) )
  else
    return ( x )
}
read.desc <- function( file, res=length(grep("res$",file)), rowskip=3 )
{
  d.idx <- if ( res ) 1 else 2
  g.idx <- if ( res ) 2 else 1

  txt <- read.delim( file, header=F, sep="\t", fill=T, skip=rowskip)[,1:2]
  desc <- as.character(txt[,d.idx])
  names(desc) <- as.character(txt[,g.idx])
  return( desc )
}
read.cls <- function( file, rowskip=2, do.lbls=(rowskip==2), verbose=F )
{
  stats <- as.vector(as.matrix(read.table(file,nrows=1,header=F)))
  if ( length(stats)!=3 )
    stop( "1st row of '.cls' file must have three (3) entries" )
  cls <- as.vector(as.matrix(read.table(file, skip=rowskip, header=F)))
  if ( length(cls)!=stats[1] )
    stop( paste("Wrong number of labels: ", length(cls), " (", stats[1], " expected)", sep="") )
  if ( my.nlevels(cls)!=stats[2] )
    stop( paste("Wrong number of label IDs: ", my.nlevels(cls), " (", stats[2], " expected)", sep="") )
  if (do.lbls)
  {
    lbls <- sapply(read.table(file,skip=1,nrows=1,comment.char=""),as.character)[-1]

    if ( length(lbls)!=stats[2] )
      stop( "label IDs must be as many as sample labels in '.cls' file" )

    # label IDs always sorted from the one associated to lowest label to highest label
    #
    levels(cls) <- lbls[order(unique(cls))]

    VERBOSE( verbose, "class labels: ", paste(levels(cls),collapse=", "), "\n" )
  }
  return( cls )    
}
########################################################################
#                         SCORES (SNR, T-SCORE)                        #
########################################################################
#
# function: FIX SD
#
fix.sd <- function( s, m, s.percent=0.2 )
{
  # function to threshold stdev as done in GeneCluster (so as to be
  # able to compare results)
  #
  if ( is.vector(s) )
  {
    if ( length(s)!=length(m) ) { stop("s and m must be same length") }

    abs.m <- abs( m )
    min.s <- abs.m * s.percent 
    min.s[min.s<s] <- s[min.s<s]
    min.s[min.s<=0] <- 0.1
    return( min.s )
  }
  else
  {
    abs.m <- abs( m )
    min.s <- s.percent * abs.m
    if ( min.s<s ) { min.s <- s }
    if ( min.s==0 ) { min.s <- 0.1 }
    return( min.s )
  }
}
# function: T.SCORE
#
t.score <- function( x, y=NULL, cls=NULL, robust=F, paired=F,
                     generalized=F, do.test=F, verbose=F, var.equal=F )
{
  # INPUT:
  #    x - m x n1 matrix (genes by experiments, condition 1)
  #    y - m x n2 matrix (genes by experiments, condition 2)
  #  OR
  #    x - m x n  matrix (genes by experiments, condition 1 & 2)
  #  cls - n vector of class labels
  #
  # OUTPUT:
  #  snr - m vector (positive are upregulated for x, or for
  #        lower label -- condition 1 -- when cls is specified)

  # some checks on the input
  #
  if ( is.null(y) & is.null(cls) )
    stop( "must specify either y or cls" )

  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2]]
    x <- x[,cls==lev[1]]
  }  
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<4 ) warning( "x has less than 4 observations\n" )
  if ( ncol(y)<4 ) warning( "y has less than 4 observations\n" )

  score <- NULL

  # paired score
  #
  if ( paired )
  {
    if ( generalized )
    {
      if ( (ncol(x)>ncol(y) & ncol(x)%%ncol(y)) |
           (ncol(y)>ncol(x) & ncol(y)%%ncol(x)) )
        stop( "x and y must be have column numbers multiple of each other\n" )

      d <- NULL

      if ( ncol(x)>ncol(y) )
        for ( i in 0:(ncol(x)/ncol(y)-1) ) {
          idx <- (i*ncol(y)+1):((i+1)*ncol(y))
          VERBOSE(verbose, "comparing x[:",idx[1],":",idx[length(idx)],
                  "] to y[1:",ncol(y),"]","\n",sep="")
          d <- cbind( d, x[,idx]-y )
        }
      else if ( ncol(y)>ncol(x) )
        for ( i in 0:(ncol(y)/ncol(x)-1) ) {
          idx <- (i*ncol(x)+1):((i+1)*ncol(x))
          VERBOSE(verbose, "comparing x[1:",ncol(x),"] to y[",
                  idx[1],":",idx[length(idx)],"]","\n",sep="")
          d <- cbind( d, x-y[,idx] )
        }
      else
        d <- x-y
    }
    else
    {
      if ( ncol(x)!=ncol(y) )
        stop( "x and y must have same number of columns\n" )

      d <- x-y
    }
    VERBOSE( verbose, "computing paired t.score .." )
    if ( robust ) {
      stop( "robust paired t.score not implemented yet" )
    }
    else {
      if (do.test) {
        score <- apply( d, 1, function(z){ Z <- t.test(z); c(Z$statistic,Z$p.value) } )
      }
      else {
        score <- apply( d, 1, function(z){ t.test(z)$statistic } )
      }
    }
    VERBOSE( verbose, " done.\n" )
    return( score )
  }
  # ELSE not paired
  #
  x.idx <- 1:ncol(x)
  TEST <- TEST <- match.fun(if(robust) wilcox.test else t.test)

  VERBOSE( verbose, ifelse( robust, "\tWilcoxon test .. ", "\tt test .. " ) )

  if(do.test) {
    score <- t(apply( cbind(x,y), 1,
                     function(z){Z <- TEST(z[x.idx],z[-x.idx],var.equal=var.equal);
                                 c(Z$statistic,Z$p.value)} ))
    colnames(score) <- c("score","p.value")

  }
  else {
    n1 <- ncol(x)
    n2 <- ncol(y)
    cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
    x <- cbind(x,y)

    s  <- x %*% cls
    s2 <- x^2 %*% cls
    s2[,1] <- (s2[,1] - (s[,1]^2)/n1) / (n1-1)
    s2[,2] <- (s2[,2] - (s[,2]^2)/n2) / (n2-1)
    s[,1] <- s[,1]/n1
    s[,2] <- s[,2]/n2

    stderr <- if (var.equal)
      sqrt( (((n1-1)*s2[,1] + (n2-1)*s2[,2])/(n1+n2-2)) * (1/n1+1/n2) )
    else
      sqrt( s2[,1]/n1 + s2[,2]/n2 )

    score <- (s[,1]-s[,2]) / stderr
  }
  VERBOSE( verbose, "done.\n" )
  return( score )
}
# function: SNR
#
snr <- function( x, y=NULL, cls=NULL, robust=F, fix=T, gc=fix,
                 paired=F, generalized=F, do.test=F, verbose=F )
{
  # INPUT:
  #    x - m x n1 matrix (genes by experiments, 1st condition)
  #    y - m x n2 matrix (genes by experiments, 2nd condition)
  #  OR
  #    x - m x n  matrix (genes by experiments, both conditions)
  #  cls - n vector of class labels
  #
  # OUTPUT:
  #  snr - m vector (positive are upregulated for x, or for
  #        lower label -- 1st condition -- when cls is specified)

  # some checks on the input
  #
  if ( is.null(y) & is.null(cls) )
    stop( "must specify either y or cls" )

  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2]]
    x <- x[,cls==lev[1]]
  }  
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<4 ) warning( "x has less than 4 observations\n" )
  if ( ncol(y)<4 ) warning( "y has less than 4 observations\n" )

  # paired SNR
  #
  if ( paired )
  {
    if ( generalized )
    {
      if ( (ncol(x)>ncol(y) & ncol(x)%%ncol(y)) |
           (ncol(y)>ncol(x) & ncol(y)%%ncol(x)) )
        stop( "x and y must have column numbers multiple of each other\n" )

      d <- NULL

      if ( ncol(x)>ncol(y) )
        for ( i in 0:(ncol(x)/ncol(y)-1) ) {
          idx <- (i*ncol(y)+1):((i+1)*ncol(y))
          VERBOSE(verbose, "comparing x[:",idx[1],":",idx[length(idx)],
                  "] to y[1:",ncol(y),"]","\n",sep="")
          d <- cbind( d, x[,idx]-y )
        }
      else if ( ncol(y)>ncol(x) )
        for ( i in 0:(ncol(y)/ncol(x)-1) ) {
          idx <- (i*ncol(x)+1):((i+1)*ncol(x))
          VERBOSE(verbose, "comparing x[1:",ncol(x),"] to y[",
                  idx[1],":",idx[length(idx)],"]","\n",sep="")
          d <- cbind( d, x-y[,idx] )
        }
      else
        d <- x-y
    }
    else
    {
      if ( ncol(x)!=ncol(y) )
        stop( "x and y must have same number of columns\n" )

      d <- x-y
    }
    VERBOSE( verbose, "computing paired SNR .." )
    if ( robust ) {
      m <- apply( d, 1, median )
      s <- apply( d, 1, mad )
    }
    else {
      m <- apply( d, 1, mean )
      s <- apply( d, 1, sd )
    }
    if ( fix ) s <- fix.sd(s,m,s.percent=.05)
    snr <- m/s
    VERBOSE( verbose, " done.\n" )

    if ( do.test )
    {
      VERBOSE( verbose, "computing asymptotic p-values (n=", ncol(d), ") ..", sep="" )
      p.value <- apply( d,1,function(z){ t.test(z)$p.value } )
      return ( snr, p.value )
      VERBOSE( verbose, " done.\n" )
    }
    return( snr )
  }
  # ELSE not paired
  #
  snr <- NULL
  if ( robust )
  {
    x.m <- apply( x, 1, median )
    x.s <- if ( gc ) apply( x, 1, sd) else apply( x, 1, mad )
    y.m <- apply( y, 1, median )
    y.s <- if ( gc ) apply( y, 1, sd ) else apply( y, 1, mad )

    # thresholding as implemented in GeneCluster (to handle zero variance)
    #
    if ( fix )
    {
      x.s <- fix.sd(x.s,x.m)
      y.s <- fix.sd(y.s,y.m)
    }
    snr <- (x.m-y.m)/(x.s+y.s)
  }
  else
  {
    n1 <- ncol(x)
    n2 <- ncol(y)
    cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
    x <- cbind(x,y)

    s  <- x %*% cls
    s2 <- x^2 %*% cls
    s2[,1] <- sqrt( (s2[,1] - (s[,1]^2)/n1) / (n1-1) )
    s2[,2] <- sqrt( (s2[,2] - (s[,2]^2)/n2) / (n2-1) )
    s[,1] <- s[,1]/n1
    s[,2] <- s[,2]/n2

    # thresholding as implemented in GeneCluster (to handle zero variance)
    #
    if ( fix ) {
      s2[,1] <- fix.sd(s2[,1],s[,1])
      s2[,2] <- fix.sd(s2[,2],s[,2])
    }    
    snr <- (s[,1]-s[,2]) / (s2[,1]+s2[,2])
  }
  if ( do.test )
  {
    VERBOSE( verbose, "computing asymptotic p-values (n=[",
             paste(ncol(x),ncol(y),sep="x"), "]) ..", sep="" )
    p.value <- apply( cbind(x,y), 1,
                      function(z){t.test(z[1:ncol(x)],z[(ncol(x)+1):length(z)],var.equal=T)$p.value} )
    return ( snr, p.value )
    VERBOSE( verbose, " done.\n" )
  }
  else {
    return( snr )
  }
}
########################################################################
#                          CLASS SHUFFLING                             #
########################################################################
#
library(survival)

many2one <- function( cls )
{
  if ( is.null(dim(cls)) || ncol(cls)==1 )
    return( as.numeric(match(drop(cls),sort(unique(drop(cls))))) )

  levs <- sort(apply(unique(cls),1,paste,collapse=""))
  cls.new <- as.numeric( match( apply(cls,1,paste,collapse=""), levs) )
  levels(cls.new) <- levs
  cls.new
}
permute.binarray.new <- function( cls, nperm=1, balanced=F, equalized=F, seed=NULL,
                                   control=NULL, verbose=F )
{
  if ( !is.null(dim(cls)) )
    if (ncol(cls)==2)
      cls <- cls[,2]
    else
      stop( "multi-dimensional class label must have 2 columns" )
  lev <- my.levels(cls)

  m  <- length(cls)
  n1 <- sum(cls==lev[1])
  n2 <- sum(cls==lev[2])
  i.max <- which.max(c(n1,n2))
  max.lev <- lev[i.max]
  min.lev <- lev[-i.max]
  idx.max <- (1:m)[cls==max.lev]
  idx.min <- (1:m)[cls==min.lev]
  n.max <- c(n1,n2)[i.max]
  n.min <- c(n1,n2)[-i.max]
  perm <- matrix( NA, nrow=nperm, ncol=length(cls), byrow=T )

  if ( !is.null(seed) ) set.seed(seed)

  if ( balanced )
  {
    perm[,idx.min] <- min.lev
    if ( length(lev)>2 )
      stop( "balanced permutation requires binary class" )
    if ( n1!=n2 && !equalized )
      stop( "balanced permutation requires same-size classes (or 'equalized=T')")
    k1 <- n.min%/%2
    k2 <- n.min%%2

    for ( i in 1:nperm )
    {
      k <- k1 + rbinom(1,1,.5)*k2
      idx.max.perm <- if ( n1!=n2 ) sample(idx.max,n.min) else idx.max
      perm[i,idx.max.perm] <- max.lev
      idx.max.perm <- sample( idx.max.perm, k )
      idx.min.perm <- sample( idx.min, 2*k-k )
      perm[i,idx.max.perm] <- min.lev
      perm[i,idx.min.perm] <- max.lev
    }
  }
  else if ( is.null(control) )
  {
    for ( i in (1:nperm) )
    {
      idx <- if (equalized) 
               c(idx.min, sample(idx.max,n.min))
             else
               c(idx.min,idx.max)

      perm[i,idx] <- sample(cls[idx],length(idx))
    }
  }
  # control for confounding (discrete) covariates
  #
  else {
    lev0 <- lev[1]
    lev1 <- lev[2]
    control <- many2one(control)
    levs <- sort(unique(control))
    idx <- lapply( levs, function(z){which(control==z)} )
    nctl <- my.tabulate( control[cls==lev1] )
    nperm.tot <- 0

    # check if enough permutations available
    #
    for ( i in 1:length(nctl) ) {
      nperm.tot <- nperm.tot + (tmp <- lchoose(length(idx[[i]]),nctl[i]))
      VERBOSE( verbose, "nchoose[",levels(control)[i],"]:   ", exp(tmp), "\n", sep="" )
    }
    VERBOSE( verbose, "nchoose[tot]:", exp(nperm.tot), "\n" )

    idx[nctl==0] <- NULL
    nctl <- nctl[nctl>0]
    if (nperm.tot<log(nperm)) 
      stop( "Number of possible permutations less than needed: ", exp(nperm.tot) )

    for ( i in 1:length(nctl) ) idx[[i]] <- c( nctl[i], idx[[i]] )

    for ( i in 1:nperm )
    {
      pidx <- unlist(sapply( idx, function(z){ sample( z[-1], z[1] ) } ))
      perm[i,] <- lev0
      perm[i,pidx] <- lev1
    }
  }
  return(perm)
}
########################################################################
#                        PERMUTATION TESTS                             #
########################################################################
#
cumineq <- function( prm, obs, dir=1, debug=F )
{
  # INPUT:
  #  - prm    n-sized array
  #  - obs    n-sized array
  # WHAT:
  #  for each entry in obs, count how many entries in prm
  #  are <= (dir=1) or >= (dir=2) than that entry
  #
  p.ord <- order(if ( dir==1 ) prm else -prm)
  o.ord <- order(if ( dir==1 ) obs else -obs)
  o.rnk <- rank(if ( dir==1 ) obs else -obs)

  # sort entries
  #
  prm <- prm[p.ord]
  obs <- obs[o.ord]

  u.obs <- unique(obs)
  cup <- c(prm,u.obs)
  cup <- cup[order(if (dir==1) cup else -cup)]
  fp <- length(cup)+1-match(obs,rev(cup))-match(obs,u.obs)

  # return values sorted according to original order (o.rnk)
  #
  return ( if (debug) 
             cbind( prm[o.rnk], obs[o.rnk], fp[o.rnk] )
           else
             fp[o.rnk] )
}
pval2fdr <- function( p )
{
  p.ord <- order(p)
  p1 <- p[p.ord]
  fdr <- p1 * length(p) / cumineq(p1,p1)
  fdr[rank(p)]
}
perm.summary <- function( obs, perm )
{
  # given an observed score vector and a permuted score vector,
  # calculate summary statistics for different types of p-values
  # this is used when "online=T" in perm.2side (to save memory)
  #
  i.p <- obs>=0
  i.n <- !i.p
  perm.srt <- sort(perm)
  #abs.ord <- order(abs(obs))

  smry <- matrix( NA, length(obs), 6 )
  smry[i.p,1] <- as.numeric( obs[i.p]<=perm[i.p] )         # 1-sided p-value
  smry[i.n,1] <- as.numeric( obs[i.n]>=perm[i.n] )         # ..
  #smry[,2]    <- as.numeric( abs(obs)<=abs(perm) )         # 2-sided p-value
  smry[,2]    <- as.numeric( obs<=perm )                   # 2-sided p-value
  smry[i.p,3] <- as.numeric( obs[i.p]<=perm.srt[i.p] )     # rank-based p-value
  smry[i.n,3] <- as.numeric( obs[i.n]>=perm.srt[i.n] )     # ..
  #smry[abs.ord,4] <-
  #  as.numeric( abs(obs)[abs.ord]<=sort(abs(perm)) )       # rank-based (2-sided) p-value
  smry[,4]    <- cumineq( abs(perm), obs=abs(obs), dir=2 ) # FPR
  smry[,5]    <- as.numeric( smry[,4]>0 )                  # FWER
  smry[,6]    <- smry[,4]                                  # FDR
  smry
}
perm.2side <- function( x, y, nperm=100, score, ngenes=NULL, seed=NULL,
                        control=NULL, balanced=F, equalized=F, match.score=F,
                        online=T, rnd=NULL, verbose=F, debug=F, ... )
{
  # Rank genes according to given score, and carry out permutation
  # test. 
  #
  # INPUT:
  #  - x           (m x n) matrix for m genes and n experiments
  #  - y           n-tuple of 'output' (can be class labels, or Surv object, etc)
  #  - nperm       number of permutation iterations
  #  - score       the score function to use to rank genes
  #  - balanced    computation based on balanced permutations (i.e., the
  #                permuted classes contain an equal proportion of the
  #                true classes)
  #  - equalized   set to true if want to carry out balanced permutation
  #                but the classes have unequal-size (not tested)
  #  - match.score if T, when balanced and unequal-sized classes, compute
  #                a matching observed score for each permuted score
  #  - online      if T, don't store permutation iterations
  #  - ...         score-specific arguments
  #
  # OUPTUT
  #  IF (online==F)
  #
  #    (m x (nperm+1)) matrix with columns
  #
  #      perm_1 perm_2 .. perm_nperm score
  #
  #    where:
  #      perm_i is the permuted scores for the i-th iteration
  #      score  is the sorted list of observed scores
  #
  # ELSE ..
  #
  #    (m x 8) matrix with columns
  #
  #      score (the observed score)
  #      p1 (1-sided p-value)
  #      etc.
  #  
  score <- match.fun( score )
  surv <- is.Surv(y)
  lev <- if ( surv ) c(0,1) else my.levels(y)
  ngenes <- min(ngenes,nrow(x))

  if ( ngenes<2 )
    stop( "ngenes must be >=2" )
  if ( length(lev)>2 )
    stop( "y must be binary" )
  if ( surv  && nrow(y)!=ncol(x) )
    stop( "y must be same length as ncol(x)" )
  if ( !surv && length(y)!=ncol(x) )
    stop( "y must be same length as ncol(x)" )

  n <- if ( is.null(dim(y)) ) length(y) else nrow(y)
  paired <- sum(y==lev[1])!=sum(y==lev[2]) && balanced && equalized && match.score

  VERBOSE( verbose, "Computing observed score .. " )
  x.obs <- score( x, cls=y, verbose=verbose, ... )
  x.idx <- c(1:round(ngenes/2), (nrow(x)-(ngenes-round(ngenes/2))+1):nrow(x))
  x.idx <- order(x.obs)[x.idx]
  x.obs <- x.obs[x.idx]
  if (any(x.obs[-1]<x.obs[-length(x.obs)])) stop("x.obs should be ordered")

  VERBOSE( verbose, "done.\n" )
  X.obs <- x.prm <- NULL
  if ( online )
  {
    x.prm <- matrix( 0, length(x.idx), 6 )
    colnames(x.prm) <- c("p1","p2","GC.p", "fpr","fwer","fdr")
  }
  else {
    x.prm <- matrix( NA, length(x.idx), nperm+1 )
    x.prm[,nperm+1] <- x.obs
    colnames(x.prm) <- c(paste( "prm", 1:nperm, sep="." ), "score" )
    X.obs <- if (paired) matrix( NA, length(x.idx), nperm )
  }
  rownames(x.prm) <- names(x.obs)

  VERBOSE( verbose, "Permutation test (", nperm, " iterations", sep="" )
  VERBOSE( balanced && verbose,  ", balanced" )
  VERBOSE( equalized && verbose, ", equalized" )
  VERBOSE( verbose, ") ..\n\t" )

  if ( !is.null(seed) )
    set.seed(seed)  
  cls.prm <- if ( surv )
    permute.binarray.new( y[,2], n=nperm, balanced=balanced, equalized=equalized,
                          control=control )
  else
    permute.binarray.new( y, n=nperm, balanced=balanced, equalized=equalized,
                          control=control )

  OBS <- x.obs
  percent <- pctstep <- max( .1, round(1/nperm,2) )
  for ( i in 1:nperm )
  {
    y.prm   <- cls.prm[i,]
    idx.prm <- !is.na(y.prm)  # needed when balanced==T and equalized==T
    y.prm   <- y.prm[idx.prm] # ..

    if ( surv ) {             # managing a survival-type labeling
      Y.prm <- y[idx.prm]
      Y.prm[y.prm==1,2] <- 1; Y.prm[y.prm==1,1] <- y[idx.prm,1][y.prm==1]
      Y.prm[y.prm==0,2] <- 0; Y.prm[y.prm==0,1] <- y[idx.prm,1][y.prm==0]
      y.prm <- Y.prm
    }
    if (paired) OBS <- score( x[,idx.prm], cls=y[idx.prm], ... )[x.idx]
    PRM <- score(x[,idx.prm], cls=y.prm, verbose=verbose, ...)[x.idx]

    if (online) {
      tmp <- perm.summary(OBS,PRM)
      if ( !all(dim(tmp)==dim(x.prm)) ) stop("incompatible dimensions for summary")
      x.prm <- x.prm + tmp
    }
    else {
      if (paired) X.obs[,i] <- OBS
      x.prm[,i] <- PRM
    }
    if ( verbose & i>=nperm*percent ) {
      VERBOSE( verbose, percent*100,"% ", sep="" )
      percent <- round(percent + pctstep,1)
    }
  }
  VERBOSE( verbose, "\n" )

  if ( online ) # complete computation of some p-values (see perm.2side.summary for 
  {             # ..the "offline" computation of these statistics
    if (length(p2.idx <- colnames(x.prm)=="p2")==0)  stop("missing 'p2' column'")
    if (length(fp.idx <- colnames(x.prm)=="fpr")==0) stop("missing 'fpr' column'")
    if (length(fd.idx <- colnames(x.prm)=="fdr")==0) stop("missing 'fdr' column'")
    x.prm <- x.prm/nperm
    x.prm[,p2.idx] <- 2*apply(cbind(x.prm[,p2.idx],1-x.prm[,p2.idx]),1,min)
    x.prm[,fp.idx] <- x.prm[,fp.idx]/nrow(x.prm)
    x.prm[,fd.idx] <- x.prm[,fd.idx]/cumineq( abs(x.obs), obs=abs(x.obs), dir=2)

    x.prm <- cbind(x.prm,fdr.test=x.prm[,fd.idx])
    x.prm[,fd.idx] <- pval2fdr(x.prm[,p2.idx])    
    x.prm <- cbind(score=x.obs,x.prm)
    if ( !is.null(rnd) )
      x.prm <- round(x.prm,rnd)
    if ( debug )
      list( x.prm=x.prm, cls.prm=cls.prm )
    else
      x.prm
  }
  else if (paired)
    list(x.prm=x.prm,x.obs=X.obs)
  else
    x.prm
}
perm.2side.summary <- function( x, rnd=NULL, probs=c(.99,.95,.50), verbose=F )
{
  # takes output of perm.2side (online=F), and compute different p-values
  #
  x.obs <- NULL
  x.avg <- NULL
  if ( is.list(x) ) # matching observed-permuted scores (when perm.2side is .. 
                    # .. called with equalized==T and match.score==T)
  {    
    if ( !all(names(x)==c("x.prm","x.obs")) )
      stop( "list(x.prm,x.obs) expected" )

    x.obs <- x$x.obs
    x.avg <- apply(x.obs,1,mean)
    x <- x$x.prm
  }
  else {
    x.obs <- matrix( x[,ncol(x)], nrow(x), ncol(x)-1 )
  }
  i.obs <- ncol(x)
  i.p <- x[,i.obs]>0
  i.n <- !i.p
  one.side <- two.side <- rnk.base <- matrix( NA, nrow(x), 1+length(probs) )
  np <- length(probs)

  # 1-sided p-value
  #
  VERBOSE( verbose, "Computing 1-sided p-value .. " )
  one.side[i.p,1] <- apply(x[i.p,-i.obs]-x.obs[i.p,]>=0, 1, sum)/(ncol(x)-1)
  one.side[i.n,1] <- apply(x[i.n,-i.obs]-x.obs[i.n,]<=0, 1, sum)/(ncol(x)-1)
  one.side[i.p,2:(1+np)] <-  t(apply(  x[i.p,-i.obs], 1, quantile, probs ))
  one.side[i.n,2:(1+np)] <- -t(apply( -x[i.n,-i.obs], 1, quantile, probs ))
  VERBOSE( verbose, "done.\n" )

  # two-sided p-value
  #
  if (F) {
  VERBOSE( verbose, "Computing 2-sided p-value .. " )
  two.side[,1]           <- apply(abs(x[,-i.obs])-abs(x.obs)>=0, 1, sum)/(ncol(x)-1)
  two.side[i.p,2:(1+np)] <-  t(apply( abs(x[i.p,-i.obs]), 1, quantile, probs ))
  two.side[i.n,2:(1+np)] <- -t(apply( abs(x[i.n,-i.obs]), 1, quantile, probs ))
  VERBOSE( verbose, "done.\n" )
  }
  if (T) {
  VERBOSE( verbose, "Computing 2-sided p-value .. " )
  two.side[,1]           <- apply(x[,-i.obs]-x.obs<=0, 1, sum)/ncol(x.obs)
  two.side[,1]           <- 2*apply( cbind(two.side[,1],1-two.side[,1]), 1, min )
  two.side[i.p,2:(1+np)] <-  t(apply( abs(x[i.p,-i.obs]), 1, quantile, probs ))
  two.side[i.n,2:(1+np)] <- -t(apply( abs(x[i.n,-i.obs]), 1, quantile, probs ))
  VERBOSE( verbose, "done.\n" )
  }
  # rank-based p-values (one-sided)
  #
  VERBOSE( verbose, "Computing rank-based p-value (one-sided) .. " )
  x[,-i.obs]   <- apply( x[,-i.obs], 2, sort )
  rnk.base[i.p,1] <- apply(x[i.p,-i.obs]-x.obs[i.p,]>=0, 1, sum)/(ncol(x)-1)
  rnk.base[i.n,1] <- apply(x[i.n,-i.obs]-x.obs[i.n,]<=0, 1, sum)/(ncol(x)-1)
  rnk.base[i.p,2:(1+np)] <-  t(apply(  x[i.p,-i.obs], 1, quantile, probs ))
  rnk.base[i.n,2:(1+np)] <- -t(apply( -x[i.n,-i.obs], 1, quantile, probs ))
  VERBOSE( verbose, "done.\n" )

  if (F) {
  # rank-based p-values (two-sided)
  #
  VERBOSE( verbose, "Computing rank-based p-value (two-sided) .. " )
  x.abs   <- apply( abs(x[,-i.obs]), 2, sort )
  abs.ord <- order( abs(x[,i.obs]))
  rnk.base2[abs.ord,1] <- apply(x.abs-abs(x.obs)[abs.ord]>=0, 1, sum )/(ncol(x)-1)
  rnk.base2[abs.ord,2:(1+np)] <- t(apply(x.abs, 1, quantile, probs))
  rnk.base2[i.n,2:(1+np)] <- -rnk.base2[i.n,2:(1+np)]
  VERBOSE( verbose, "done.\n" )
  }
  # empirical FPR
  #
  VERBOSE( verbose, "Computing FPR .. " )
  fps <- matrix( NA, nrow(x), ncol(x.obs) )
  for ( i in 1:ncol(x.obs) ) {
    fps[,i] <- cumineq( abs(x[,i]), obs=abs(x.obs[,i]), dir=2 )
  }
  fpr <- apply( fps, 1, mean )/nrow(x)
  VERBOSE( verbose, "done.\n" )

  # empirical FWER
  #
  VERBOSE( verbose, "Computing FWER .. " )
  fwer <- apply( fps>0, 1, sum )/ncol(fps)
  VERBOSE( verbose, "done.\n" )

  # empirical FDR
  #
  VERBOSE( verbose, "Computing FDER .. " )
  fdr <-  ( fpr*length(fpr) ) / cumineq( abs(x[,i.obs]), obs=abs(x[,i.obs]),dir=2)
  VERBOSE( verbose, "done.\n" )

  colnames(one.side) <- colnames(two.side) <- colnames(rnk.base) <- 
    c( "p.val", paste(round(probs*100,1),"%",sep="") )

  rownames(one.side) <-
    rownames(two.side) <-
      rownames(rnk.base) <-
          names(fpr) <-
            names(fwer) <-
              names(fdr) <- rownames(x)
  out <- list(obs=x[,i.obs], one.side=one.side, two.side=two.side,
              rnk.base=rnk.base, fpr=fpr, fwer=fwer, fdr=fdr, avg=x.avg )
  if ( !is.null(rnd) )
    out <- lapply(out, round, rnd )
  perm.2side.pvalues(out)
}
perm.2side.pvalues <- function( obj, rnd=Inf )
{
  if ( is.null(obj$obs) )         stop( "obs missing from obj" )
  if ( is.null(obj$one.side) )    stop( "one.side missing from obj" )
  if ( is.null(obj$two.side) )    stop( "two.side missing from obj" )
  if ( is.null(obj$rnk.base) )    stop( "rnk.base missing from obj" )
  #if ( is.null(obj$rnk.twoside) ) stop( "rnk.twoside missing from obj" )
  if ( is.null(obj$fpr) )         stop( "fpr missing from obj" )
  if ( is.null(obj$fwer) )        stop( "fwer missing from obj" )
  if ( is.null(obj$fdr) )         stop( "fdr missing from obj" )

  round(cbind(score=obj$obs, p1=obj$one.side[,1], p2=obj$two.side[,1],
              GC.p=obj$rnk.base[,1], fpr=obj$fpr, fwer=obj$fwer, fdr=obj$fdr),rnd)
}
