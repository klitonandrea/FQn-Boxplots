library("pcaPP")
library("MASS")
library("ICSNP")
library("aplpack")
library("graphics")
library("ggplot2")
library("robcor")

M_1_SQRT2 <- 1 / sqrt(2)

#' Correlation function
#'
#' Calculate correlation coefficient through FQn estimates
#' @param x a vector 
#' @param y a vector
#' @return correlation coefficient and scale estimates
#' @export
#' @example 
#' CorrFQn(rnorm(50), rnorm(50))
CorrFQn <- function(x, y)
{
  # robcor(x, y)
  FQnX <- fqn(x)
  FQnY <- fqn(y)
  
  u <- (x/FQnX + y/FQnY) * M_1_SQRT2
  v <- (x/FQnX - y/FQnY) * M_1_SQRT2
  
  FQnU <- fqn(u)
  FQnV <- fqn(v)
  
  FQnU <- FQnU^2
  FQnV <- FQnV^2
  
  corr <- (FQnU - FQnV)/(FQnU + FQnV)
  
  return(list(Corr = corr, VarX = FQnX, VarY = FQnY))
}

#' Location estimate function
#'
#' Estimate location using different estimates
#' @param x a vector
#' @param y a vector
#' @return pair of values (c1, c2)
#' @export
#' @example 
#' EstimateLocation(rnorm(100), rnorm(100))
EstimateLocation <- function(x, y, estimatorType="CompMedian"){
  require(pcaPP)
  require(ICSNP)
  require(depth)
  
  if ( estimatorType == "CompMedian" )
    c(median(x), median(y))
	else
    if ( estimatorType == "SpatialMedian" )
      spatial.median(matrix(c(x,y), ncol = 2))#  HR.Mest(matrix(c(x,y), ncol=2))
	else
    if ( estimatorType == "HR")
    {
      val <- HR.Mest(matrix(c(x,y), ncol = 2))
      return(val$center)
    }
    else
      if ( estimatorType == "Tukey")
      {
        depthMedian <- try(med(matrix(c(ns$X, ns$Y), ncol = 2), method = "Tukey") )
        if ( (depthMedian == "try-error")[1] )
          depthMedian <- c(NA, NA)
        else
          depthMedian <- depthMedian$median
        return(depthMedian)
      }
  else
	    c(mean(x), mean(y))
}

#' Scale estimate function
#' 
#' Estimate scale in 2 dimensions
#' @param x vector
#' @param y vector
#' @param estimatorType 
#' @return pair of values for x and y respectively
EstimateScale <- function(x, y, estimatorType="Qn"){
  #if ( estimatorType == "Qn")
    c(fqn(x), fqn(y))
}

#' Distance function
#' 
#' Estimates distances between 2 points given their covariation matrix
#' @param from a point in multivariate dimension
#' @param to a point in multivariate dimension
#' @param CovMatr covariation matrix
#' @return scalar, the distance between points
#' @export
EstimateDistance <- function(from, to, CovMatr, estimatorType="Mahalanobis"){
  # TODO
  if ( estimatorType == "Mahalanobis")
    {
      sqrt(mahalanobis(from, to, CovMatr))
    ##mahalanobis(from, to, CovMat$cov)
    }
    else
      sqrt(mahalanobis(from, to, diag(ncol(from))))
  }

#' Create covariation matrix
#' 
#' Create covariation matrix given the datapoints and their center
#' @param x vector of scalars
#' @param y vector of scalars
#' @param center 2-dimensional center
#' @export
#' @return covariation matrix
EstimateCovMatr <- function(x, y, center)
{
  rroFQn <- CorrFQn(x - center[1], y - center[2])
  c2     <- rroFQn$Corr * rroFQn$VarX * rroFQn$VarY
  
  retVal <- matrix(c(rroFQn$VarX^2, c2, c2, rroFQn$VarY^2), ncol = 2)
  
  return(retVal)
}

#' Calculate hinge
#' 
#' @param x the first dimension of the points in 2-dimensional space
#' @param y the second dimension of the point in 2-dimensional space
#' @return half of the points closer to the center
GetHalfPointsFromCenterLocation <- function(x, y, center){
  ##center <- rep(center, times=length(x))
  
  dists <- EstimateDistance(matrix(c(x, y), ncol = 2), center, EstimateCovMatr(x, y, center))
  sortres <- sort(dists, index.return = TRUE)
  #fix(`sortres`)
  sortres$ix[1:(length(sortres$ix) %/% 2)]
  }
  
#' Calculate fence points
#' 
#' Calculate fence points - returns logic vector
#' @param x
#' @param y
#' @param center
#' @param scaleX
#' @param scaleY
#' @param alpha parameter alpha for fence, default 1.0
#' @param  whichInds vector containing the indices for the points 
#' @return logical vector indicating which points are in the fence
GetFencePoints <- function(x, y, center, scaleX, scaleY, alphaF=1.0, whichInds, type="ellipse"){
  ##dists <- EstimateDistance(c(x, y), center)
  inds    <- !(1:length(x) %in% whichInds)
  xmin    <- min(x[inds]); xmax <- max(x[inds])
  ymin    <- min(y[inds]); ymax <- max(y[inds])
  factor  <- alphaF
  left    <- max(xmin, (center[1] - factor*scaleX))
  bottom  <- max(ymin, (center[2] - factor*scaleY))
  right   <- min(xmax, (center[1] + factor*scaleX))
  up      <- min(ymax, (center[2] + factor*scaleY))
  ncenter <- c((right + left)/2.0, (up + bottom)/2.0)
  #print(c(left, bottom, right, up, factor))
  if ( type == "ellipse") {
    newCoordX <- x - ncenter[1]
    newCoordY <- y - ncenter[2]
    R   <- sqrt(newCoordX^2 + newCoordY^2)
    phi <- atan2(newCoordY, newCoordX)
    
    A <- (right - left)/2.0  #factor*scaleX
    B <- (up - bottom)/2.0  #factor*scaleY
    
    RE <- A * B/sqrt((B*cos(phi))^2 + (A*sin(phi))^2)
    rV <- ((R - RE)/RE <= 0.001)
    
    return( rV & inds)
  }
  else
    (abs(x - center[1]) < 5*scaleX*alphaF) & (abs(y - center[2]) < 5*scaleY*alphaF)
}
  
#' Calculate FQn boxplot
#' 
#' main function to calculate bi-variate boxplot
#' @param x a numeric vector or matrix
#' @param y optional numeric vector
#' @param alpha the fence coefficient
#' @return a list with boxplot calculated parameters
#' @export
FQnBoxplot <- function(x, y, alpha=1.5, verbose = FALSE){
  require(MASS)
  ## converting vectors x, y to a matrix
  dataMatrix <- NA
  if ( is.matrix(x) ) {
    dataMatrix <- x
  }else{
    dataMatrix <- matrix(c(x, y), length(x), 2)
  }

  if ( verbose == TRUE)
  {
    plot(dataMatrix)
    View(dataMatrix)
  }
    
  fqnCorr <- CorrFQn(dataMatrix[,1], dataMatrix[,2])
  
  sigmaX <- fqnCorr$VarX#CovMat$cov[1,1]
  sigmaY <- fqnCorr$VarY#CovMat$cov[2,2]
  
  ##calculate angle of rotation
  phi <- atan(2*fqnCorr$Corr*sigmaX*sigmaY / (sigmaX^2 - sigmaY^2)) / 2
  #print(phi)
  ##build rotation matrix
  rotationMatrix <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2)
  boxplotCenter <- EstimateLocation(dataMatrix[,1], dataMatrix[,2], estimatorType = "SpatialMedian")
  ## rotate axis
  #View(dataMatrix)
  dataprim <-  rotationMatrix %*% t(matrix(c(dataMatrix[,1] - boxplotCenter[1], dataMatrix[,2] - boxplotCenter[2]), ncol = 2 ) )
  dataprim <- t(dataprim)
	#if(verbose == TRUE)
    #plot(c(0,0), xlab="X'", ylab="Y'", xlim=c(-6, 6), ylim=c(-15,15))
  xprim <- dataprim[,1]
  yprim <- dataprim[,2]

  if ( verbose == TRUE )
    plot(dataprim, xlab = "X'", ylab = "Y'")
  
  ##calculate center and scale of rotated data
  boxplotCenter <- EstimateLocation(xprim, yprim, estimatorType = "SpatialMedian")
  boxplotScale <- EstimateScale(xprim, yprim) ##c(sqrt(sigmaX), sqrt(sigmaY)) ##
  
  ## find 50% of points closer to center
  closerpoints <- GetHalfPointsFromCenterLocation(xprim, yprim, boxplotCenter)
  ##View(closerpoints)
  ## build convex hull
	xfence <- xprim[closerpoints]
	yfence <- yprim[closerpoints]
  hinge <- chull(xfence, yfence)
  hinge <- matrix(c(c(xfence[hinge], xfence[hinge[1]]), c(yfence[hinge], yfence[hinge[1]])), ncol = 2)
  
  #DX <- max(c(abs(min(hinge[,1] - boxplotCenter[1])), abs(max(hinge[,1] - boxplotCenter[1]))))
  #DY <- max(c(abs(min(hinge[,2] - boxplotCenter[2])), abs(max(hinge[,2] - boxplotCenter[2]))))
  
   
  ## build fence
	fenceInd <- GetFencePoints(xprim, yprim, boxplotCenter, boxplotScale[1], boxplotScale[2], alphaF = alpha, closerpoints)
  	fencepoints <- dataprim[c(which(fenceInd), closerpoints),] #matrix(c(xprim[fenceInd], yprim[fenceInd]), ncol=2)
  
	fence <- chull(fencepoints)
  
	fencepoints <- rbind(fencepoints[fence,], fencepoints[fence[1],])

  ## get outliers
	outliersInds <- which(!(1:nrow(dataMatrix) %in% c(which(fenceInd), closerpoints)))
	outliers <- dataMatrix[outliersInds,]
  
  ## rotate back the points
	rotationMatrix <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), ncol = 2)
	boxplotCenter <- t(rotationMatrix %*% t(matrix(boxplotCenter, 1, 2)))
	hinge <- t(rotationMatrix %*% t(hinge))
	fence <- t(rotationMatrix %*% t(fencepoints))
	bpret <- list(Correlation  = fqnCorr, 
	              Phi          = phi, 
	              Points       = dataMatrix, 
	              Center       = boxplotCenter, 
	              BoxplotScale = boxplotScale, 
	              Hinge        = hinge, 
	              Fence        = fence, 
	              Outlier      = outliers, 
	              OutlierInds  = outliersInds)
  class(bpret) <- "FQnBoxplot"
  return(bpret)
  }

#' Plot a list of bivariate FQn boxplots
#' 
#' @param bps list of boxplots
#' @param names boxplot naming, by default the boxplots will be enumerated
#' @param xLabel the x label
#' @param yLabel the y label
#' @export
plotBPs <- function(bps, names = c(), xLabel = "", yLabel = ""){
  library("ggplot2")
  N <- length(bps)
  if ( N == 0)
    return
  
  if ( length(names) == 0)
    names <- 1:N
  
  categories <- factor(c("Outlier","Fence","Hinge", "Center", "Points"))
  
  bpData <- data.frame(
    BoxplotId     = c(),
    Category      = c(),
    CategoryValue = c(),
    X             = c(),
    Y             = c(), 
    BPName        = c()
    )
  
  #names(bpData)<-c("BoxplotId",
  #  "Category",
  #  "CategoryValue",
  #  "X",
  #  "Y"
  #  )
  for (i in 1:N) {
    #add Outliers
    outlength <- length(bps[[i]]$Outlier)
    if (outlength > 0)
    {
      outm <- bps[[i]]$Outlier
      
      if (outlength == 2)
        outm <- matrix(bps[[i]]$Outlier, ncol = 2)
      
    outs <- data.frame(
      BoxplotId     = rep(names[i], outlength),
      Category      = rep(categories[1], outlength),
      CategoryValue = rep(i + 0.1, outlength),
      X             = outm[,1],
      Y             = outm[,2],
      BPName        = names[i]
      )
    bpData <- rbind(bpData, outs)
  }
    #add Fence
    fencelength <- nrow(bps[[i]]$Fence)
    print(fencelength)
    fences <- data.frame(
      BoxplotId     = rep(names[i], fencelength),
      Category      = rep(categories[2], fencelength),
      CategoryValue = rep(i + 0.2, fencelength),
      X             = bps[[i]]$Fence[,1],
      Y             = bps[[i]]$Fence[,2],
      BPName        = names[i]
      )
    bpData <- rbind(bpData, fences)
    #add Hinge
    hingelength <- nrow(bps[[i]]$Hinge)
    hinges <- data.frame(
      BoxplotId     = rep(names[i], hingelength),
      Category      = rep(categories[3], hingelength),
      CategoryValue = rep(i + 0.3, hingelength),
      X             = bps[[i]]$Hinge[,1],
      Y             = bps[[i]]$Hinge[,2],
      BPName        = names[i]
      )
    bpData <- rbind(bpData, hinges)
    #add Center
    center <- data.frame(
      BoxplotId     = names[i],
      Category      = categories[4],
      CategoryValue = i + 0.4,
      X             = bps[[i]]$Center[1],
      Y             = bps[[i]]$Center[2],
      BPName        = names[i]
      )
    #Merge last boxplot
    bpData <- rbind(bpData, center)
  }
  
  aD <- subset(bpData, Category == "Fence")
  aD$BoxplotId = rep("ALL", length(aD$BoxplotId));
  map <- c(geom_polygon(aes(x = X, y = Y, group = CategoryValue, fill = BPName, alpha = 0.5), data = aD))
  
#View(bpData)
  (p <- ggplot(bpData, aes(x = X, y = Y)) 
   + geom_polygon(data = subset(bpData, Category == "Fence"), fill = '#A52A2A') #
   + geom_polygon(data = subset(bpData, Category == "Hinge"), fill = 'green', alpha = 0.5)#alpha=0.5 aes(fill=factor(Category) )
   #+ geom_point(data=subset(bpData, Category=="Outlier"), colour = 'black') #aes(fill=factor(Category) )
   + geom_point(data = subset(bpData, Category == "Center"), colour = 'black'))#aes(fill=factor(Category) )
  
   if (length(which(bpData$Category == "Outlier")) > 0)
     p <- (p + geom_point(data = subset(bpData, Category == "Outlier"), colour = 'red'))
  if (xLabel != "")
    p <- p + xlab(xLabel)
  if (yLabel !=  "")
    p <- p + ylab(yLabel)
   (p + facet_wrap(~ BoxplotId, nrow = 1))  + map + scale_fill_hue(name = "Categories", breaks = names) + scale_alpha_continuous(guide = FALSE) #+ opts(legend.position="none") scales="free_x",
}

#' Plot a bivariate FQn-boxplot
#' 
#' @param bp boxplot
#' @param show.outliers by default set to true
#' @export
plot2dBP <- function(bp, show.outliers=TRUE, ordernum=1)
{
 
  #CovariationMCD = CovMat, Phi=phi, Points=dataMatrix, Center=boxplotCenter, 
  #BoxplotScale=boxplotScale, 
  #Hinge=hinge, Fence = fence, Outlier=outliers, OutlierInds=outliersInds
  
  
  positionsFence <- data.frame(
    x = bp$Fence[,1],
    y = bp$Fence[,2]
    )
  #datapoly<-merge(rep(5, times=nrows(bp$Fence)), )
  (p <- ggplot(positionsFence, aes(x = x, y = y)) + geom_polygon(aes(fill = 5 )))
  
  positionsHinge <- data.frame(
    x = bp$Hinge[,1],
    y = bp$Hinge[,2]
    )
  p <- p + geom_polygon(data = positionsHinge, fill = 5.1, alpha = 0.5)
  if (show.outliers == TRUE) {
    outliersPos <- data.frame(x = bp$Outlier[,1], y = bp$Outlier[,2])
    p           <- p + geom_point(data = outliersPos)
  }
  p + geom_point(data = data.frame(x = bp$Center[1], y = bp$Center[2]), colour = 5.2) +
    theme(legend.position = "none")
}

#' Another plotting method
#' 
#' @export
PlotBoxplot <- function(box, scatterPlot=TRUE, add=FALSE)
{
	#x1 <- box$Points[,1]
	#x2 <- box$Points[,2]
  
  ## plot outliers
  if (!add)
    plot(box$Points[,1], box$Points[,2], xlab = "X", ylab = "Y", bty = "n", type = "n") 
  
  polygon(box$Fence, border = 'gray82', col = 'gray82')
  if (is.array(box$Outlier)) {
    points(box$Outlier, pch = "*")
    ##plot fence
  #lines(box$Fence)
  }
   # else
    #{
     ##}
	##plot hinge
  #a<-seq(from=0, to=2*pi, by=pi/6)
  #ex<-box$BoxplotScale[1]*cos(a-box$Phi)+box$Center[1]
  #ey<-box$BoxplotScale[2]*sin(a-box$Phi)+box$Center[2]
  #lines(ex, ey)
    polygon(box$Hinge, border = 'gray57', col = 'gray57')
	#lines(box$Hinge)
	
	
	if (scatterPlot)
		points(box$Points)
	##plot center
	points(box$Center, pch = "o")
  
	
}



Bagplot <- function(sig, alpha = 2.7){ 
#   signalX<-sig$Signal[, 1]
#     noiseX<-sig$Noise[,1] #if(nrow(sig$Noise)>1) sig$Noise[2:nrow(sig$Noise), 1] else c()
#     signalY<-sig$Signal[ , 2]
#     noiseY<- sig$Noise[, 2]
#     x<-c(signalX, noiseX)
#     y<-c(signalY, noiseY)
  dt <- sig
  if (is.matrix(sig)) {
    dt <- list(X = sig[, 1], Y = sig[, 2])
  }
  oldSeed <- .Random.seed; 
    bg <- compute.bagplot(dt$X, dt$Y, factor = alpha, approx.limit = 100000)
  # Push random seed to old state 
      assign(".Random.seed", oldSeed, envir = globalenv())
  return(bg)
}

TestSet <- function()
{
  t <- data.frame(
    X = c(-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 0, 0, 0, 0, 0, 0, 0, -0.25, -0.25, 0.25, 0.25),
    Y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, -0.5, -0.25, -0.22, 0, 0.22, 0.25, 0.5, 0.25, -0.25, 0.25, -0.25)
  )
  
  return(t)
}

