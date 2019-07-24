# basic plotting functions for assessing output

require('fields')

temp.colors=function(n=25){
      m <- floor(n/2)
        blues <- hsv(h=.65, s=seq(1,0,length=m+1)[1:m])
        reds <- hsv(h=0, s=seq(1,0,length=m+1)[1:m])
        c(blues,if(n%%2!=0) "#FFFFFF", reds[m:1])
    }


f.clr=function(vec){
  return(log(vec)-mean(log(vec)))
}

# rIts: space, taxon, time, iteration?
# rMean: space, taxon, time (mean across iterations? )
# p: taxon number?
spatialEOFbyTaxon=function(p,rMean,rIts,transform='asin'){
  # this returns EOFs for a specific taxon where the patterns are with respect to the locations (so true eofs) and there is a score for each time for each of the nT patterns (where nT is the number of time points)
  # since we work with a single taxon, the clr transformation doesn't make sense as the transformation is based on having all P taxa, and if one applied it and then took only the transformed data for the taxon of interet, one will have different data values even if the abundance of beech, say, were 24% in two different time periods because the abundunace of the other taxa changed between the two times
  # scoreMean is the mean of the scores on the patterns (eofs) for each time point, averaging over the scores for each iteration; each row is one pattern and each column a time point
  # scoreUpper is the 97.5th percentile and scoreLower the 2.5th percentile of the scores for a given time and pattern; each row is one pattern and each column a time point
  # scoreIts is the entire set of scores by iteration: patterns X times X iterations
  # scoresOnMeans are the scores calculated on rMean rather than by iteration - note that the EOFs are computed based on rMean so that the patterns are the same for every iteration; otherwise the meaning of the scores would change between iterations; each row is one pattern and each column a time point
  # loadings are the weights on the locations for each pattern - each column is a pattern and each row is a grid cell
  # eigenvals are the variances associated with each pattern with higher variance meaning more variability explained by that pattern
  nT=dim(rIts)[3]
  S=dim(rIts)[1]
  nIts=dim(rIts)[4]

  calcScores=function(data){
    return((lm(data~mysvd$u-1))$coef)
  }

 
  data=asin(sqrt(rMean[,p,]))

  locMeans=rowMeans(data)
  data=data-locMeans 

  mysvd=svd(data)

  workItsMat=matrix(c(asin(sqrt(rIts[,p,,]))),S,nT*nIts)-locMeans

  tmp=lm(workItsMat~mysvd$u-1)$coef
  scoreIts=array(tmp,c(nT,nT,nIts))  # there are only 31 patterns because only 31 time points, so this is the score on each pattern for each time point for each iteration
  scoreMean=apply(scoreIts,c(1,2),mean)  # mean score on each pattern for each time point
  scoreUpper=apply(scoreIts,c(1,2),quantile,.975)  # quantiles of scores on each pattern for each time point
  scoreLower=apply(scoreIts,c(1,2),quantile,.025)
  
  scoresOnMeans=diag(mysvd$d)%*%t(mysvd$v)
  
  return(list(scoreMean=scoreMean,scoreUpper=scoreUpper,scoreLower=scoreLower,scoreIts=scoreIts,loadings=mysvd$u,eigenvals=mysvd$d,scoresOnMeans=scoresOnMeans))
}

spatialEOFallTaxa=function(rMean,rIts,transform='clr'){ # could also do transform='asin' or transform='none'
  # this returns spatial EOFs where they are calculated based on all taxa
  # the patterns are with respect to the locations (so true eofs) and there is a score for taxon for each time for each of the S=256 patterns
  # one can either invoke the clr transformation that deals with having proportions or one could invoke 'asin' for the arcsin sqrt transformation. I assume that the EOFs based on the covariance are desired, but to change this, just change the line below that does the svd (mysvd=svd(....
  # scoreMean is the mean of the scores on the patterns (eofs) for each taxon and time point, averaging over the scores for each iteration; each row is one pattern and each column a time point
  # scoreUpper is the 97.5th percentile and scoreLower the 2.5th percentile of the scores for a given taxon and time; each row is one pattern
  # scoreIts is the entire set of scores by iteration: patterns X taxa X times X iterations
  # scoresOnMeans are the scores calculated on rMean rather than by iteration - note that the EOFs are computed based on rMean so that the patterns are the same for every iteration; otherwise the meaning of the scores would change between iterations; each row is one pattern
  # loadings are the weights on the locations for each pattern - each column is a pattern and each row is a grid cell
  # eigenvals are the variances associated with each pattern with higher variance meaning more variability explained by that pattern
  nT   = dim(rIts)[3]
  S    = dim(rIts)[1]
  P    = dim(rIts)[2]
  nIts = dim(rIts)[4]

  calcScores=function(data){
    return((lm(data~mysvd$u-1))$coef)
  }
  
  print('Got dimensions!')
  
  if(transform=='asin'){
    data=asin(sqrt(rMean))
    workItsMat=matrix(c(asin(sqrt(rIts))),S,P*nT*nIts)
  }
  
  if(transform=='clr'){
    tmp=apply(rMean,c(1,3),f.clr)
    data=array(NA,c(S,P,nT))
    for(p in 1:P){
      data[,p,]=tmp[p,,]
    }
    workIts=apply(rIts,c(1,3,4),f.clr)
    tmp=array(NA,c(S,P,nT,nIts))
    for(p in 1:P){
      tmp[,p,,]=workIts[p,,,]
    }
    workItsMat=matrix(c(tmp),S,P*nT*nIts)
  }
  
  print('Transformed!')
  
  tmp=data
  data=matrix(c(data),S,P*nT)

  locMeans=rowMeans(data)
  data=data-locMeans 

  mysvd=svd(data)

  workItsMat=workItsMat-locMeans

  tmp        = lm(workItsMat~mysvd$u-1)$coef
  scoreIts   = array(tmp,c(S,P,nT,nIts))  # there are 256 patterns, and this is the score on each pattern for each taxon and each time point for each iteration
  scoreMean  = apply(scoreIts,c(1,2,3),mean)  # mean score on each pattern for each taxon and time point
  scoreUpper = apply(scoreIts,c(1,2,3),quantile,.975)  # quantiles of scores on each pattern for each taxon and time point
  scoreLower = apply(scoreIts,c(1,2,3),quantile,.025)
  
  scoresOnMeans=array(diag(mysvd$d)%*%t(mysvd$v),c(S,P,nT))
  
  return(list(scoreMean=scoreMean,scoreUpper=scoreUpper,scoreLower=scoreLower,scoreIts=scoreIts,loadings=mysvd$u,eigenvals=mysvd$d,scoresOnMeans=scoresOnMeans))
}


taxaEOFs=function(rMean,rIts){
  # this returns EOFs where the patterns are with respect to the taxa and there is a score for each grid cell at each time for each of the 9 patterns (the 10th is meaningless because of sum to one constraint
  # scoreMean is the mean of the scores on the patterns (eofs) for each grid cell by time, averaging over the scores for each iteration
  # scoreUpper is the 97.5th percentile and scoreLower the 2.5th percentile of the scores for a given cell and time and pattern
  # scoreIts is the entire set of scores by iteration
  # scoresOnMeans are the scores calculated on rMean rather than by iteration - note that the EOFs are computed based on rMean so that the patterns are the same for every iteration; otherwise the meaning of the scores would change between iterations
  # loadings are the weights on the taxa for each pattern - each column is a pattern and each row is a taxon
  # eigenvals are the variances associated with each pattern with higher variance meaning more variability explained by that pattern
  nT=dim(rIts)[3]
  S=dim(rIts)[1]
  nIts=dim(rIts)[4]
  P=dim(rIts)[2]
  mymat=matrix(0,nr=S*nT,nc=P)
  for(p in 1:P){
    mymat[,p]=c(rMean[,p,])
  }
  
#   f.clr=function(vec){
#     return(log(vec)-mean(log(vec)))
#   }
  

  data=apply(mymat,1,f.clr) # transformation to account for proportional data
  taxaMeans=rowMeans(data)
  data=data-taxaMeans # this puts rMean output into a P=10 by n= S*nT matrix; n is the number of 'observations' or 'replicates'

  mysvd=svd(data)

  calcScores=function(data){
    return((lm(data~mysvd$u-1))$coef)
  }

  workIts=array(NA,c(P,S,nT,nIts))
  for(p in 1:P){
    workIts[p,,,]=rIts[,p,,]
  }
  workItsMat=matrix(c(workIts),P,S*nT*nIts)
  workItsMat=apply(workItsMat,2,f.clr)-taxaMeans

  tmp=lm(workItsMat~mysvd$u-1)$coef

  scoreIts=array(tmp,c(P,S,nT,nIts))
  scoreMean=apply(scoreIts,c(1,2,3),mean)
  scoreUpper=apply(scoreIts,c(1,2,3),quantile,.975)
  scoreLower=apply(scoreIts,c(1,2,3),quantile,.025)
  
  scoresOnMeans=array(diag(mysvd$d)%*%t(mysvd$v),c(P,S,nT))  # scores on the rMean values at each location and time

  tmp1=tmp2=tmp3=tmp4=array(NA,c(S,P,nT))
  tmp5=array(NA,c(S,P,nT,nIts))
  for(p in 1:P){
    tmp1[,p,]=scoreMean[p,,]
    tmp2[,p,]=scoreLower[p,,]
    tmp3[,p,]=scoreUpper[p,,]
    tmp4[,p,]=scoresOnMeans[p,,]
    tmp5[,p,,]=scoreIts[p,,,]
  }
  
  return(list(scoreMean=tmp1,scoreUpper=tmp3,scoreLower=tmp2,scoreIts=tmp5,loadings=mysvd$u,eigenvals=mysvd$d,scoresOnMeans=tmp4))
}


addPondLocs=function(tNew,tOld=NULL){
# plots black circles for ponds with pollen data and green open circles for ponds without
# use tOld argument if you want to plot information for two time points
  points(pondLocs[pondsUsed[[tNew]],],pch=16,cex=.8)
  points(pondLocs[pondsNotUsed[[tNew]],],col='green')
  if(!is.null(tOld)){
    points(pondLocs[pondsUsed[[tOld]],],pch=16,cex=.8)
    points(pondLocs[pondsNotUsed[[tOld]],],col='green')
  }
}


taxaSignif=function(p1,p2,rIts,locs,restrict=FALSE,inputRestricted=FALSE,ylab='',xlab=''){
# compares significance of differences in composition between two predictions
# red colors indicate more of the taxon in the older time period; blue in the newer
  nComparisons=length(rIts[1,1,])
  less=rIts[,p1,]<rIts[,p2,]
  less=matrix(apply(less,1,sum),sqrt(S),sqrt(S))
  more=nComparisons-less
  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }
  
  xs=sort(unique(locs[,1]))
  ys=sort(unique(locs[,2]))[rows]
  image(xs,ys,more[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,col=tim.colors(32)[21:32],xaxt='n',yaxt='n')
  image(xs,ys,less[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,col=tim.colors(32)[12:1],add=TRUE,xaxt='n',yaxt='n')
  
  lines(ma$x,ma$y,lwd=1)
  lines(ct$x,ct$y,lwd=1)
  lines(ri$x,ri$y,lwd=1)
  lines(ny$x,ny$y,lwd=1)
  lines(nh$x,nh$y,lwd=1)
}

timeSignifSynch=function(p1,p2,rItsNew,rItsOld,locs,restrict=FALSE,inputRestricted=FALSE,ylab='',xlab=''){
# compares synchrony of significance of differences in composition between two predictions in time and two taxa
# blue colors indicate both taxa show lower abundance in the newer time period; red both higher abundance in the newer; green that the first taxa shows higher in the older and second higher in the newer, and brown that the first shows higher in the new and second higher in the older
  nComparisons=length(rItsNew[1,1,])

  bothLess=rItsNew[,p1,]<rItsOld[,p1,] & rItsNew[,p2,]<rItsOld[,p2,]
  bothLess=matrix(apply(bothLess,1,sum),sqrt(S),sqrt(S))
  bothMore=rItsNew[,p1,]>rItsOld[,p1,] & rItsNew[,p2,]>rItsOld[,p2,]
  bothMore=matrix(apply(bothMore,1,sum),sqrt(S),sqrt(S))
  firstLess=rItsNew[,p1,]<rItsOld[,p1,] & rItsNew[,p2,]>rItsOld[,p2,]
  firstLess=matrix(apply(firstLess,1,sum),sqrt(S),sqrt(S))
  secondLess=rItsNew[,p1,]>rItsOld[,p1,] & rItsNew[,p2,]<rItsOld[,p2,]
  secondLess=matrix(apply(secondLess,1,sum),sqrt(S),sqrt(S))

  #more=nComparisons-less

  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }
  
  xs=sort(unique(locs[,1]))
  ys=sort(unique(locs[,2]))[rows]
  more=nComparisons-bothLess
  image(xs,ys,bothMore[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=tim.colors(32)[21:32])
  image(xs,ys,bothLess[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=tim.colors(32)[12:1],add=TRUE)

  ramp <- colorRamp(c("lightgreen", "darkgreen"))
  greens=rgb( ramp(seq(0, 1, length =12)), max = 255)

  ramp <- colorRamp(c("wheat", "chocolate4"))
  browns=rgb( ramp(seq(0, 1, length =12)), max = 255)

  image(xs,ys,firstLess[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,col=greens,add=T)
  image(xs,ys,secondLess[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,col=browns,add=TRUE)
  
  lines(ma$x,ma$y,lwd=1)
  lines(ct$x,ct$y,lwd=1)
  lines(ri$x,ri$y,lwd=1)
  lines(ny$x,ny$y,lwd=1)
  lines(nh$x,nh$y,lwd=1)
}


timeDiffs=function(p,rItsNew,rItsOld,locs,lowCut=0.9,highCut=0.95,maxDiff=0.225,nColors=4,restrict=FALSE,inputRestricted=FALSE,ylab='',xlab='',nIncr=50,lwd=2,legend=TRUE){
# compares differences in composition between two predictions, at two levels of significance
# red colors indicate more of the taxon in the newer time period; red in the newer
# hashing indicates lower level of significance
# it's good to make maxDiff=(nColors*2+1)*value/2 where value is a roundish number such as .05,.10 that is the size of the interval for each color in the color bar, e.g., maxDiff=.225 when nColors=4 and you want intervals of .05, which gives the intervals (-.225,-.175),(-.175,-.125),(-.125,-.075),(-.075,-.025),(-.025,.025),...,(.175,.225)

  # nColors*2+1 is the number of different colors plotted
  # maxDiff is the maximum difference in abundance plotted; other values are thresholded to maxDiff and negative maxDiff
  # nIncr and lwd control the diagonal hashing (poorly...)
  # lowCut and highCut are the posterior probability cutoffs
  nComparisons=length(rItsNew[1,1,])
  less=rItsNew[,p,]<rItsOld[,p,]
  less=matrix(apply(less,1,sum),sqrt(S),sqrt(S))
  more=nComparisons-less
  meanNew=rowMeans(rItsNew[,p,])
  meanOld=rowMeans(rItsOld[,p,])
  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }

  lessProp=less/nComparisons
  moreProp=more/nComparisons
  lowSignif=lessProp>=lowCut | moreProp>=lowCut
  highSignif=lessProp>=highCut | moreProp>=highCut

  diff=meanNew-meanOld
  diff[!lowSignif]=0
  diff=matrix(thresh(diff,maxDiff,-maxDiff),sqrt(S),sqrt(S))
  
  xs=sort(unique(locs[,1]))
  ys=sort(unique(locs[,2]))[rows]
  if(legend){
    image.plot(xs,ys,diff[,rows],zlim=c(-maxDiff,maxDiff),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=temp.colors(nColors*2+1))
  } else{
    image(xs,ys,diff[,rows],zlim=c(-maxDiff,maxDiff),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=temp.colors(nColors*2+1))
  }
  xincr=diff(range(xs))/(2*length(xs))
  yincr=diff(range(ys))/(2*length(ys))
  xrg=range(xs)+c(-1,1)*xincr
  yrg=range(ys)+c(-1,1)*yincr

  incr=(yrg[2]-xrg[1]-yrg[1]+xrg[2])/nIncr
  for(i in 1:nIncr){
    abline(a=yrg[1]-xrg[2]+incr*(i-1),b=1,col='white',lwd=2,lty=1)
  }
  box()
  diff[!highSignif]=NA
  if(legend){
    largeplot=image.plot.cp(xs,ys,diff[,rows],zlim=c(-maxDiff,maxDiff),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=temp.colors(nColors*2+1),add=T)
    par(plt=largeplot)
  } else{
    image(xs,ys,diff[,rows],zlim=c(-maxDiff,maxDiff),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=temp.colors(nColors*2+1),add=T)
  }
  lines(ma$x,ma$y,lwd=1)
  lines(ct$x,ct$y,lwd=1)
  lines(ri$x,ri$y,lwd=1)
  lines(ny$x,ny$y,lwd=1)
  lines(nh$x,nh$y,lwd=1)
  box()
}

image.plot.cp=function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 0.05, graphics.reset = FALSE, horizontal = FALSE, 
    offset = 2 * legend.width, bigplot = NULL, smallplot = NULL, 
    legend.only = FALSE, col = tim.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) 
        big.plot <- old.par$plt
    if (legend.only) 
        graphics.reset <- TRUE
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, horizontal = horizontal, 
        offset = offset, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        image(..., add = add, col = col)
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    temp <- list(...)
    iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    ix <- 1
    if (!horizontal) {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
        axis(4, mgp = c(3, 1, 0), las = 2)
        box()
    }
    else {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
            col = col)
        box()
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        return(bigplot)
        #invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        return(bigplot)
        #invisible()
    }
}


timeSignif=function(p,rItsNew,rItsOld,locs,restrict=FALSE,inputRestricted=FALSE,ylab='',xlab=''){
# compares significance of differences in composition between two predictions
# red colors indicate more of the taxon in the newer time period; blue in the older
  nComparisons=length(rItsNew[1,1,])
  less=rItsNew[,p,]<rItsOld[,p,]
  less=matrix(apply(less,1,sum),sqrt(S),sqrt(S))
  more=nComparisons-less
  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }
   
  xs=sort(unique(locs[,1]))
  ys=sort(unique(locs[,2]))[rows]
  image(xs,ys,more[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=tim.colors(32)[21:32])
  image(xs,ys,less[,rows],zlim=c(.9*nComparisons,nComparisons),xlab=xlab,ylab=ylab,xaxt='n',yaxt='n',col=tim.colors(32)[12:1],add=TRUE)
  
  lines(ma$x,ma$y,lwd=1)
  lines(ct$x,ct$y,lwd=1)
  lines(ri$x,ri$y,lwd=1)
  lines(ny$x,ny$y,lwd=1)
  lines(nh$x,nh$y,lwd=1)
}

feature=function(p,less,locs,restrict=FALSE,inputRestricted=FALSE,nComparisons=10000){
  # plots feature significance for one prediction for one taxon
  par.old=par()
  less=less[[p]]
  more=nComparisons-less
  diag(more)=0
  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }

  par(mfrow=c(m2,m1),mai=c(0,0,0,0),xaxt='n',yaxt='n')
  for(row in rowseq){
    for(col in 1:16){
      image(1:m1,1:m2,matrix(more[((row-1)*16+col),][used],m1,m2),zlim=c(.9*nComparisons,nComparisons),xlab='',ylab='',col=tim.colors(32)[21:32],cex.axis=2,yaxt='n',xaxt='n')
      image(1:m1,1:m2,matrix(less[((row-1)*16+col),][used],m1,m2),zlim=c(.9*nComparisons,nComparisons),xlab='',ylab='',col=tim.colors(32)[12:1],cex.axis=2,yaxt='n',xaxt='n',add=TRUE)
      text(col,row-min(rowseq)+1.001,'x',cex=1.4)  # can't figure out why this is not printing 'x's on bottom row
    }
  }
  par(par.old)
}



# modify this to deal better with input format
timePlot=function(times,mean,lower=NULL,upper=NULL,ylim=NULL,xlab='',ylab=''){
  # plots time series with confidence bars (e.g., to plot coefficient values over time)
  par(mfrow=c(2,5))
  for(p in 1:ncol(mean)){
    plot(-1*times,mean[,p],type='l',ylim=range(c(lower,upper,mean)),ylab=ylab,xlab=xlab)
    points(-1*times,mean[,p],pch=16,cex=.8)
    lines(-1*times,lower[,p],col=2,lty=2)
    lines(-1*times,upper[,p],col=2,lty=2)
    abline(h=0,col='gray',lty=3)
    title(taxa[p])
  }
}

timePlotDecomp=function(times,mean,lower=NULL,upper=NULL,its=NULL,ylim=NULL,xlab='',ylab=''){
  # plots time series with confidence bars (e.g., to plot coefficient values over time)
  mns=its
  tmp=apply(its,c(2,3),mean)
  for(t in 1:length(times)){
    mns[t,,]=tmp
  }
  its=its-mns
  devMean=apply(its,c(1,2),mean)
  devLower=apply(its,c(1,2),quantile,.025)
  devUpper=apply(its,c(1,2),quantile,.975)

  par(mfrow=c(2,5))
  for(p in 1:ncol(mean)){
    plot(-1*times,mean[,p],type='n',ylim=range(c(lower,upper,mean)),ylab=ylab,xlab=xlab)
#    points(-1*times,mean[,p],pch=16,cex=.8)
    polygon(c(rev(-1*times),-1*times),c(rev(lower[,p]),upper[,p]),col='gray80',border='gray80')
    lines(-1*times,mean[,p])
    lines(-1*times,devLower[,p]+mean(mean[,p]),col=2)
    lines(-1*times,devUpper[,p]+mean(mean[,p]),col=2)
    abline(h=0,col='gray50',lty=3)
    title(taxa[p])
  }
}

# modify to deal better with input format
vegDiagram=function(cell,times){
  # analogous to a pollen diagram - plots vegetation predictions over time
  mean=matrix(0,nr=length(times),nc=P)
  lower=upper=mean
  for(tt in 1:length(times)){
    mean[tt,]=eval(as.name(paste('rMean',times[tt],sep='')))[cell,]
    lower[tt,]=eval(as.name(paste('rLow',times[tt],sep='')))[cell,]
    upper[tt,]=eval(as.name(paste('rUp',times[tt],sep='')))[cell,]
  }
  maxvals=apply(upper,2,max)
  maxvals=signif(maxvals,digits=1)+.1
  for(p in 1:P){
    if(maxvals[p]<.2){maxvals[p]=.2}
  }
  par(mfrow=c(1,P),mai=c(0.4,0.15,0.4,0.05),mgp=c(0,0.3,0))
  p=1
  plot(mean[,p],-1*times,type='p',xlim=c(0,maxvals[p]),pch=16,cex=.3,xlab='',ylab='',tck=-.04)
  title(taxaNames[p])
  for(vv in seq(.1,.9,by=.1)){
    abline(v=vv,lty=2,col='grey')
  }
  lines(mean[,p],-1*times,pch=16,cex=0.3)
  lines(lower[,p],-1*times,col=2)
  lines(upper[,p],-1*times,col=2)
  par(mai=c(0.4,0.15,0.4,0.05))
  for(p in 2:P){
    plot(mean[,p],-1*times,type='p',xlim=c(0,maxvals[p]),yaxt='n',pch=16,cex=.3,xlab='',ylab='',tck=-.04)
    title(taxaNames[p])
    for(vv in seq(.1,.9,by=.1)){
      abline(v=vv,lty=2,col='grey')
    }
    lines(mean[,p],-1*times)
    lines(lower[,p],-1*times,col=2)
    lines(upper[,p],-1*times,col=2)
  }
}
  
vegPlusPollenDiagram=function(cell,times,pollenAges,pollenProp,lwd=1){
  # analogous to a pollen diagram - plots vegetation predictions over time with pollen overlaid; pollen might go outside figure limits...
  # pollenProp should have rows equal ages and columns be the 10 taxa with the data being proportions
  mean=matrix(0,nr=length(times),nc=P)
  lower=upper=mean
  for(tt in 1:length(times)){
    mean[tt,]=eval(as.name(paste('rMean',times[tt],sep='')))[cell,]
    lower[tt,]=eval(as.name(paste('rLow',times[tt],sep='')))[cell,]
    upper[tt,]=eval(as.name(paste('rUp',times[tt],sep='')))[cell,]
  }
  maxvals=apply(upper,2,max)
  maxvals=signif(maxvals,digits=1)+.1
  for(p in 1:P){
    if(maxvals[p]<.2){maxvals[p]=.2}
  }
  par(mfrow=c(1,P),mai=c(0.4,0.15,0.4,0.05),mgp=c(0,0.3,0))
  p=1
  plot(mean[,p],-1*times,type='p',xlim=c(0,maxvals[p]),pch=16,cex=.3,xlab='',ylab='',tck=-.04)
  title(taxaNames[p])
  for(vv in seq(.1,.9,by=.1)){
    abline(v=vv,lty=2,col='grey',lwd=lwd)
  }
  lines(mean[,p],-1*times,pch=16,cex=0.3,lwd=lwd)
  lines(lower[,p],-1*times,col=2,lwd=lwd)
  lines(upper[,p],-1*times,col=2,lwd=lwd)
  lines(pollenProp[,p],pollenAges,col='blue',lwd=lwd)
  par(mai=c(0.4,0.15,0.4,0.05))
  for(p in 2:P){
    plot(mean[,p],-1*times,type='p',xlim=c(0,maxvals[p]),yaxt='n',pch=16,cex=.3,xlab='',ylab='',tck=-.04)
    title(taxaNames[p])
    for(vv in seq(.1,.9,by=.1)){
      abline(v=vv,lty=2,col='grey',lwd=lwd)
    }
    lines(mean[,p],-1*times,lwd=lwd)
    lines(lower[,p],-1*times,col=2,lwd=lwd)
    lines(upper[,p],-1*times,col=2,lwd=lwd)
    lines(pollenProp[,p],pollenAges,col='blue',lwd=lwd)
  }
}

vegPlusPollenDiagramDecomp=function(cell,timeInd,times,pollenAges,pollenProp,maxvals=NULL,lwd=1){
  # analogous to a pollen diagram - plots vegetation predictions over time with pollen overlaid but now with long-term avg + deviations; pollen might go outside figure limits...
  # pollenProp should have rows equal ages and columns be the 10 taxa with the data being proportions

  smps=rIts[cell,,timeInd,]
  lterm=matrix(0,nr=P,nc=3)
  lterm[,1]=apply(smps,1,mean)
  lterm[,2]=apply(smps,1,quantile,c(.025))
  lterm[,3]=apply(smps,1,quantile,c(.975))

  mns=smps
  for(t in 1:length(times)){
    mns[,t,]=apply(smps,c(1,3),mean)
  }
  smps=smps-mns
  tv=array(0,c(P,length(times),3))
  tv[,,1]=apply(smps,c(1,2),mean)
  tv[,,2]=apply(smps,c(1,2),quantile,.025)
  tv[,,3]=apply(smps,c(1,2),quantile,.975)

  mean=rMean[cell,,timeInd]
  lower=rLow[cell,,timeInd]
  upper=rUp[cell,,timeInd]

  if(is.null(maxvals)){
    pollMax=apply(tmp,2,max)
    maxvals=apply(upper,1,max)
    maxvals=apply(cbind(pollMax,maxvals),1,max)
    maxvals=signif(maxvals,digits=1)+.1
    for(p in 1:P){
      if(maxvals[p]<.2){maxvals[p]=.2}
    }
  }
  par(mfrow=c(1,P),mai=c(0.4,0.15,0.4,0.05),mgp=c(0,0.3,0))
  p=1
  plot(mean[p,],-1*times,type='n',xlim=c(0,maxvals[p]),pch=16,cex=.3,xlab='',ylab='',tck=-.04)
  title(taxaNames[p])
#  lines(lower[,p],-1*times,col=2)
#  lines(upper[,p],-1*times,col=2)
  polygon(c(rev(lower[p,]),upper[p,]),c(rev(-1*times),-1*times),col='gray80',border='gray80')
  for(vv in seq(.1,.9,by=.1)){
    abline(v=vv,lty=2,col='gray50',lwd=lwd)
  }
  lines(mean[p,],-1*times,pch=16,cex=0.3,lwd=lwd)
#  points(mean[p,],-1*times,pch=16,cex=0.3)
  lines(lterm[p,1]+tv[p,,1],-1*times,pch=16,cex=0.3,lwd=lwd)
  lines(lterm[p,1]+tv[p,,2],-1*times,col=2,lwd=lwd)
  lines(lterm[p,1]+tv[p,,3],-1*times,col=2,lwd=lwd)
  lines(pollenProp[,p],pollenAges,col='blue',lwd=lwd)
  points(pollenProp[,p],pollenAges,col='blue',pch=16,cex=.3)
  par(mai=c(0.4,0.15,0.4,0.05))
  for(p in 2:P){
    plot(mean[p,],-1*times,type='n',xlim=c(0,maxvals[p]),yaxt='n',pch=16,cex=.3,xlab='',ylab='',tck=-.04)
    title(taxaNames[p])
    polygon(c(rev(lower[p,]),upper[p,]),c(rev(-1*times),-1*times),col='gray80',border='gray80')
     for(vv in seq(.1,.9,by=.1)){
      abline(v=vv,lty=2,col='grey50',lwd=lwd)
    }
    lines(mean[p,],-1*times,lwd=lwd)
   #  points(mean[p,],-1*times,pch=16,cex=0.3)
   #lines(lower[,p],-1*times,col=2)
    #lines(upper[,p],-1*times,col=2)
    lines(lterm[p,1]+tv[p,,1],-1*times,lwd=lwd)
    lines(lterm[p,1]+tv[p,,2],-1*times,col=2,lwd=lwd)
    lines(lterm[p,1]+tv[p,,3],-1*times,col=2,lwd=lwd)
    lines(pollenProp[,p],pollenAges,col='blue',lwd=lwd)
    points(pollenProp[,p],pollenAges,col='blue',pch=16,cex=.3)
  }

}

thresh=function(vec,up=NULL,lo=NULL){
  if(!is.null(up)){
    vec[vec>up]=up
  }
  if(!is.null(lo)){
    vec[vec<lo]=lo
  }
  return(vec)
}

surfMap=function(p,rMat,locs,max=0.4,min=0,restrict=FALSE,inputRestricted=FALSE,legend=TRUE,statelwd=2){
  # surface map for a single taxon
  if(!restrict){
    used=1:S
    m1=m2=sqrt(S)
    rows=1:16
  } else{
    if(inputRestricted){
      m1=16;m2=9;
      used=1:(m1*m2)
      rows=1:9
    } else{
      used=49:(256-64)
      m1=16;m2=9
      rows=4:12
    }
  }

  tmp=list()
  tmp$x=unique(locs[used,1])
  tmp$y=unique(locs[used,2])
  m1=length(tmp$x)
  m2=length(tmp$y)
  tmp$z=matrix(rMat[used,p],m1,m2)
  tmp$z=thresh(tmp$z,max)
  if(legend){
    image.plot(tmp$x,tmp$y,tmp$z,xaxt='n',yaxt='n',col=tim.colors(32),zlim=c(min,max),xlab="",ylab="")
  } else{
    image(tmp$x,tmp$y,tmp$z,col=tim.colors(32),xaxt='n',yaxt='n',xlab="",ylab="",zlim=c(min,max))
  }
  lines(ma$x,ma$y,lwd=statelwd)
  lines(ct$x,ct$y,lwd=statelwd)
  lines(ri$x,ri$y,lwd=statelwd)
  lines(nh$x,nh$y,lwd=statelwd)
  lines(ny$x,ny$y,lwd=statelwd)
}


pieMap=function(proportions,centers,restrict=FALSE,inputRestricted=FALSE,xlim=c(629000-6000,809000+6000),ylim=c(4621000-6000,4801000+6000),radius=NULL,scale=1,xlab='x',ylab='y',...){
# plots multiple pie composition charts as a map
  if(!restrict){
    used=1:nrow(centers)
  } else{ # assumes subset of grid
    if(inputRestricted){
      used=1:144
    } else{
      used=49:(256-64)
    }
  }
  centers=as.matrix(centers[used,])
  proportions=as.matrix(proportions[used,])
  if(is.null(xlim)){
    rg=range(centers[,1])
    df=(scale-1)*diff(range(centers[,1]))
    xlim=c(rg[1]-df,rg[2]+df)
  }
  if(is.null(ylim)){
    rg=range(centers[,2])
    df=(scale-1)*diff(range(centers[,2]))
    ylim=c(rg[1]-df,rg[2]+df)
  }
  plot(centers,type='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  n=length(centers[,1])
  if(is.null(radius)){
    radius=.025*diff(range(centers[,1]))
  }
  minVal=min(proportions)
 # eps=1e-10*max(proportions)
 # if(minVal==0){
 #   warning("Some proportions are zero; proceeding with jittering.")
 #   proportions=proportions+eps
 # }
  if(minVal<0){
    stop("Some proportions are less than zero.")
  }
  if(length(radius)==1){ radius=rep(radius,n)}
  for(i in 1:n){
    if(sum(proportions[i,])>0){
      minVal=min(proportions[i,])
      if(minVal==0){
        warning("Some proportions are zero; proceeding with jittering.")
        eps=1e-10*max(proportions[i,])
        proportions[i,]=proportions[i,]+eps
      }
      pieAdd(as.vector(proportions[i,]),as.vector(centers[i,]),radius=radius[i],...)
    } else{
      points(centers[i,],pch='X')
    }
  }
  lines(ma$x,ma$y,lwd=1)
  lines(ct$x,ct$y,lwd=1)
  lines(ri$x,ri$y,lwd=1)
  lines(ny$x,ny$y,lwd=1)
  lines(nh$x,nh$y,lwd=1)
}


pieAdd= function (x, center, labels = names(x), edges = 200, radius = 0.8, density = NULL, angle = 45, col = NULL, border = NULL, lty = NULL, ...) # modified from the pie() function in R
{
    if (!is.numeric(x) || any(is.na(x) | x <= 0)) 
        stop("pie: `x' values must be positive.")
    if (is.null(labels)) 
        labels <- rep("",length(x))
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)

    pin <- par("pin")
    nx <- length(dx)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "black","lightblue", "red","darkblue","yellow",
              "purple","orange","lightgreen","darkgreen")
        else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    for (i in 1:nx) {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
        xc <- c(cos(t2p), 0) * radius + center[1]
        yc <- c(sin(t2p), 0) * radius + center[2]
        polygon(xc, yc, density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i],...)
        t2p <- 2 * pi * mean(x[i + 0:1])
        xc <- cos(t2p) * radius + center[1]
        yc <- sin(t2p) * radius + center[2]
        if (!is.na(lab <- labels[i]) && lab != "") {
            lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
            text(1.1 * xc, 1.1 * yc, lab, xpd = TRUE, adj = ifelse(xc < 
                0, 1, 0), ...)
        }
    }
    invisible(NULL)
}


# template code for using above functions to analyze output

if(FALSE){

# piemap of predictions for third time point (500 in 300-3000 predictions), restricting to approximate region with ponds
pieMap(rMean[,,3],gridLocs,restrict=TRUE)
# piemap of aggregated pollen for third time point
pieMap(cMat[,,3],pondLocs)

# surface maps for third time point
surfMap(1,rMean[,,3],gridLocs,max=0.6,restrict=TRUE)

par(mfrow=c(3,3),mar=c(2.1,2.1,3.1,1.1)) 
  for(p in 1:(P-1)){
    surfMap(p,rMean[,,3],gridLocs,max=0.6,restrict=TRUE)
    title(taxa[p])
  }

# feature significance map

feature(2,less[[3]],gridLocs,restrict=TRUE,nComparisons=nComparisons) # 2nd taxon, 3rd time period
# ignore the 'graphical parameter' warning messages



# maps of significant differences between two time periods (the third and fourth time periods) for 2nd taxon

# pure significance
timeSignif(2,rIts[,,3,],rIts[,,4,],gridLocs,restrict=TRUE)
# actual mean differences with only significant ones shown
timeDiffs(p=1,rIts[,,3,],rIts[,,4,],gridLocs,restrict=T,nColors=4,maxDiff=.225)

# putting all taxa on a page
 par(mfrow=c(3,3),mar=c(2.1,2.1,3.1,1.1)) 
  for(p in 1:(P-1)){
    timeSignif(p,rIts[,,3,],rIts[,,4,],gridLocs,restrict=TRUE)
    title(taxa[p])
  }

# map of significant sychrony between two taxa (pine and oak) between two time periods (the third and fourth time periods)
timeSignifSynch(8,7,rIts[,,3,],rIts[,,4,],gridLocs,restrict=TRUE)




# time series plots of coefficients
whichts=1:23 # use 300-2500 time interval (if the output is for 300-3000 time period)
timePlotDecomp(times[whichts],b1sMean[whichts,],b1sLow[whichts,],b1sUp[whichts,],b1s[whichts,,])

# maps of significant differences between two taxa (2nd and 7th here) for 3rd time period
taxaSignif(2,7,rIts[,,3,],gridLocs,restrict=TRUE)


# vegetation + pollen diagam for 300-2500 time interval for cell 139 (Snake Pond I think)
dat=read.table('PollenTimeSeries.csv',sep=',',header=T)
i=20
subdat=dat[dat$sitenumber==i,]
age=-subdat$'cal.age'
tmp=subdat[,5:14]
tmp=tmp/apply(tmp,1,sum)
whichTimeIndices=1:23
vegPlusPollenDiagramDecomp(139,whichTimeIndices,seq(300,2500,by=100),age,tmp)

# example calculation of taxa eofs and uncertainty

myTaxaEOFs=taxaEOFs(rMean,rItsSub)
# some plotting ideas
# map of scores on first EOF for 12th time point
surfMap(1,myTaxaEOFs$scoreMean[,,12],gridLocs,restrict=T)
# time series plot of first EOF at 80th grid point for all time points, with uncertainty
tsplot(myTaxaEOFs$scoreMean[80,1,],ylim=c(-12,12),xlab='time',ylab='scores')
lines(myTaxaEOFs$scoreUpper[80,1,],lty=2)
lines(myTaxaEOFs$scoreLower[80,1,],lty=2)

# example of spatial EOFs specific to a taxon

beechEOF=spatialEOFbyTaxon(1,rMean,rItsSub)
# plot the first eof
surfMap(1,beechEOF$loadings,gridLocs,max=0.2,min=-0.2,restrict=T)    # here '1' is the pattern desired not the taxon
# plot time series and uncertainty of scores for first eof
tsplot(beechEOF$scoreMean[1,],ylim=c(-4,4),xlab='time',ylab='scores')
lines(beechEOF$scoreUpper[1,],lty=2)
lines(beechEOF$scoreLower[1,],lty=2)

# spatial EOFs for all taxa
allEOF=spatialEOFallTaxa(rMean,rItsSub,transform='clr')
# plot the first eof
surfMap(1,allEOF$loadings,gridLocs,max=0.2,min=-0.2,restrict=T)    # here '1' is the pattern desired not the taxon
# plot time series and uncertainty of scores for first eof, second taxon (birch)
tsplot(allEOF$scoreMean[1,2,],ylim=c(-60,60),xlab='time',ylab='scores')
lines(allEOF$scoreUpper[1,2,],lty=2)
lines(allEOF$scoreLower[1,2,],lty=2)
# plot time series of posterior mean scores for first eof, all taxa
for(p in 1:P){
  if(p==1){
    tsplot(allEOF$scoreMean[1,p,],ylim=c(-40,40),xlab='time',ylab='scores')
  } else{
    lines(allEOF$scoreMean[1,p,],col=p)
  }
}

}
