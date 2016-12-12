
MonteCarloLocation <- function(nsize = 10000)
{
  library("depth")
  noiseType <- sample(1:11, nsize, replace=TRUE)
  locX <- rep(0, nsize)
  locY <- rep(0, nsize)
  
  sampleSize <- sample(c(50, 100, 1000, 500, 200, 300, 10000), nsize, replace=TRUE)
  
  parTable <- data.frame(
    NT = noiseType,
    SS = sampleSize
    )
  
  resultTable <- data.frame(
    LocX = locX,
    LocY = locY,
    NoiseType = noiseType,
    SampleSize = sampleSize,
    Median = rep(0, nsize),
    Spatial = rep(0, nsize),
    HR = rep(0, nsize),
    DepthTukey = rep(0, nsize),
    DepthLiu = rep(0, nsize),
    DepthOja = rep(0, nsize),
    DepthSpat = rep(0, nsize)
    )
 
  processF <- function(x){
    #print(c(x["NT"], x["SS"]))
    ns<-GenerateNoiseByModel(type=x["NT"], N=100) #N=x["SS"])
    if( length(ns$X) != 100 )
      print("data size not equal 100")
    medianEst<-EstimateLocation(ns$X, ns$Y, estimatorType = "CompMedian")
    medianE <- medianEst[1]^2 + medianEst[2]^2
    
    spatmedianEst <- NA
    HREst <- NA
    HREstE <- NA
    spatmedianEst<-try(EstimateLocation(ns$X, ns$Y, estimatorType = "SpatialMedian"))
    if( (spatmedianEst == "try-error")[1] )
      spatmedianEst <- c(NA, NA)
    spatmedianE <- spatmedianEst[1]^2 + spatmedianEst[2]^2
    
    HREst<-try(EstimateLocation(ns$X, ns$Y, estimatorType = "HR") )
    if( (HREst == "try-error")[1] )
      HREst <- c(NA, NA)
    
    if(length(HREst) != 2)
      print(HREst)
    else
      HREstE <- HREst[1]^2 + HREst[2]^2
    
    depthMedian <- try(med(matrix(c(ns$X, ns$Y), ncol = 2), method="Tukey") )
    if( (depthMedian == "try-error")[1] )
      depthMedian <- c(NA, NA)
    else
      depthMedian <- depthMedian$median
    depthMedianE <- depthMedian[1]^2 + depthMedian[2]^2
    
    depthMedianL <- try(med(matrix(c(ns$X, ns$Y), ncol = 2), method="Liu") )
    if( (depthMedianL == "try-error")[1] )
      depthMedianL <- c(NA, NA)
    else
      depthMedianL <- depthMedianL$median
    depthMedianLE <- depthMedianL[1]^2 + depthMedianL[2]^2
    
    depthMedianO <- try(med(matrix(c(ns$X, ns$Y), ncol = 2), method="Oja") )
    if((depthMedianO == "try-error")[1] )
      depthMedianO <- c(NA, NA)
    else
      depthMedianO <- depthMedianO$median
    depthMedianOE <- depthMedianO[1]^2 + depthMedianO[2]^2
      
    depthMedianS <- try(med(matrix(c(ns$X, ns$Y), ncol = 2), method="Spatial") )
    if((depthMedianS == "try-error")[1] )
      depthMedianS <- c(NA, NA)
    else
      depthMedianS <- depthMedianS$median
    depthMedianSE <- depthMedianS[1]^2 + depthMedianS[2]^2
    #print(c(medianE, spatmedianE, HREstE, depthMedianE, depthMedianLE, depthMedianOE, depthMedianSE))
    c(medianE, spatmedianE, HREstE, depthMedianE, depthMedianLE, depthMedianOE, depthMedianSE)
  }
  
  result<-apply(parTable, 1, processF)
  #resultTable$Median <- sqrt( result[1,])
  
  #resultTable$Spatial <- sqrt( result[3,]^2 + result[4,]^2)
  #resultTable$HR <- sqrt( result[5,]^2 + result[6,]^2)
  
  return(result)
}