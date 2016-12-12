EqualPoints<-function(x1, y1, x2, y2){
    x1==x2 & y1==y2
}

GetOutliers<-function(bgplot, data){
  bgOut<-try(bgplot$pxy.outlier)
  if(class(bgplot)!="bagplot")
  {
    fix(bgplot)
  }
  Out<-ifelse(length(data$Noise)>0, 
    data$Noise[!is.na(data$Noise[, 1]) & !is.na(data$Noise[, 2]), ],
    matrix(, 0, 2) )
  
  if(class(Out)=="numeric")
    Out<-matrix(Out, ncol=2)
  TN<-0
  FN<-0
  if((length(bgOut)>0) && (nrow(Out)>0)){
    for(i in 1:nrow(Out)){
      for(j in 1:nrow(bgOut)){
        if(EqualPoints(Out[i,1], Out[i,2], bgOut[j,1], bgOut[j,2]))
          TN<-TN+1
        
      }
    }
  }
  
  nbgOut<-0
  nOut<-0
  if(!is.null(nrow(bgOut)))
    nbgOut<-nrow(bgOut)
    
  if(!is.null(nrow(Out))) nOut<-nrow(Out)
  FN<- nbgOut - TN  
  FP<-nOut - TN 
  TP<-nrow(data$Signal)-FN
  sensitivity<-0
  if((TP+FN)>0)
    sensitivity<-TP/(TP+FN)
  specificity<-0
  if((FP+TN)>0)
    specificity<-TN/(FP+TN)
  param<-c(TN, FN, FP, TP, sensitivity, specificity)
  return(param)
}

ROCBagplot<-function(data){
  alpha <- 1
  outlierCount<-1
  counter<-1
  ROC<-matrix(c(0,0,0), 1, 3)
  
  while(outlierCount>0)
  {
    bg<-Bagplot(alpha, data)
    param<-GetOutliers(bg, data)
    outlierCount<-param[1]
    ##fix(param)
    ROC<-rbind(ROC, matrix(c(alpha, param[5], param[6]), 1, 3))
    if(param[5]==1.0) break
    alpha<-alpha+0.1
    counter<-counter+1
    if(alpha > 100){ 
      print(c(alpha, param[5], param[6]))
      break
    }
  }
  
  ROC

}

SimulateMonteCarloBPOld<-function(num){
  Simulation<-list()
  
  for(i in 1:num){
      Simulation[i]<-list(ROCBagplot())
    }
    Simulation
}

GetConfusionMatrix<-function(SignalNoise, a=2.7, method="elliptic"){
  signalX<-SignalNoise$Signal[, 1]
  noiseX<-if(length(SignalNoise$Noise)>0) SignalNoise$Noise[!is.na(SignalNoise$Noise[, 1]), 1] else c()
  signalY<-SignalNoise$Signal[, 2]
  noiseY<- if(length(SignalNoise$Noise)>0) SignalNoise$Noise[!is.na(SignalNoise$Noise[, 2]), 2] else c()
  x<-SignalNoise$X
	y<-SignalNoise$Y
  SigThresh<-nrow(SignalNoise$Signal) 
  ROC<-matrix(c(0,0,0), 1, 3)
  #print(SigThresh)
  
  boxPlot <-c()
  if(method == "elliptic")
	  boxPlot <- BuildBoxplot(x, y, alpha=a)
  else
    if(method == "mve")
      boxPlot <- mvemethod(matrix(c(x, y), ncol = 2 ))
  else
    if(method == "project")
      boxPlot <- promethod(data.frame(X = x, Y = y))
  else
    if(method == "mgv")
      boxPlot <- mgvmethod(matrix(c( x, y), ncol = 2, byrow=TRUE) )
  else
    if(method == "mcd")
      boxPlot <- mcdmethod(data.frame(X = x, Y = y))
  
  #PlotBoxplot(boxPlot, scatterPlot=FALSE)
  outlierCount<-length(boxPlot$OutlierInds)
  TN<-which(boxPlot$OutlierInds>SigThresh)
  FN<-which(boxPlot$OutlierInds<=SigThresh)
  #View(boxPlot$OutlierInds)
  TN<-length(TN)
  
  FN<-length(FN)
  FP<- length(noiseX) - TN
  TP<-length(signalX)-FN
  sensitivity<-0
  if((TP+FN)>0)
  sensitivity<-TP/(TP+FN)
  specificity<-0
  if((FP+TN)>0)
    specificity<-TN/(FP+TN)
  else
    specificity <- 1
  #print(c(TN, FN, FP, TP))
  param<-c(TN, FN, FP, TP, sensitivity, specificity)
  return(param)
}

GetROCPoints<-function(SignalNoise){
  signalX<-SignalNoise$Signal[, 1]
  noiseX<-SignalNoise$Noise[, 1]
  signalY<-SignalNoise$Signal[, 2]
  noiseY<- SignalNoise$Noise[, 2] 
	x<-c(signalX, noiseX)
	y<-c(signalY, noiseY)
  SigThresh<-nrow(SignalNoise$Signal) - 1
  alpha <- 0.1
  outlierCount<-1
  counter<-1
  ROC<-matrix(c(0,0,0), 1, 3)
  
  while(outlierCount>0)
  {
	  boxPlot<-BuildBoxplot(x, y, alpha)
    outlierCount<-length(boxPlot$OutlierInds)
    TN<-which(boxPlot$OutlierInds>SigThresh)
    FN<-which(boxPlot$OutlierInds<=SigThresh)
    TN<-length(TN);
    FN<-length(FN);
    FP<-nrow(SignalNoise$Noise) - TN-1
    TP<-nrow(SignalNoise$Signal)-FN-1
    sensitivity<-TP/(TP+FN)
    specificity<-TN/(FP+TN)
    param<-c(TN, FN, FP, TP, sensitivity, specificity)
    ##fix(param)
    ROC<-rbind(ROC, matrix(c(alpha, sensitivity, specificity), 1, 3))
    if(sensitivity==1.0) break
    alpha<-alpha+0.1
    counter<-counter+1
  }
  
  ROC
}

SimulateMonteCarlo<-function(num){
  Simulation<-list()
  result<-list()
  counter<-1
  AUCs<-c()
  for(j in 10:10)
  {
    
    for(i in 1:num){
      print(c(i, j))
        data<-GenerateNoiseByModel(type=j)
        rx<-try(GetROCPoints(data))
        if(class(GetROCPoints)=="try-error")
        {
          fix(j)
          fix(data)
        }
        AUCs[i]<-AUC(rx)
        Simulation[i]<-list(rx)
      }
      ##result[counter]<-list(try(GetMedianForEachAlpha(Simulation)))
      result[counter]<-list(AUCs)
      counter<-counter+1
  }
  result
}

SimulateMonteCarloBP<-function(num){
  Simulation<-list()
  result<-list()
  counter<-1
  AUCs<-c()
  for(j in 10:10)
  {
    d<-list()
    for(i in 1:num){
      
        d[i]<-list(GenerateNoiseByModel(type=j))
    }
    ##fix(d)
    for(i in 1:num){
        print(c(i, j))
        rx<-try(ROCBagplot(d[[i]]))
        if(class(ROCBagplot)=="try-error")
        {
          fix(j)
          fix(data)
        }
        AUCs[i]<-AUC(rx)
        Simulation[i]<-list(rx)
        print(rx)
      }
      ##result[counter]<-list(try(GetMedianForEachAlpha(Simulation)))
      result[counter]<-list(AUCs)
      counter<-counter + 1
  }
  result
}


GetMedianForEachAlpha<-function(rocs){
  N<-length(rocs)
  tempMax<-rocs[[1]][,1]
  for(i in 2:N){
    tempMax<-pmax(tempMax, rocs[[i]][,1])
    }
  MaxLength<-length(tempMax)
  
  data<-matrix(c(0, 0,0), 1, 3)
  #fix(tempMax)
  for(i in 2:MaxLength){
    p<-matrix(c(0, 0), 1, 2)
    for(j in 1:N){
      if(i<= length(rocs[[j]][,1])){
        if(!is.na(rocs[[j]][i,2]) && !is.na(rocs[[j]][i,3]))
        {
        pp<-matrix(c(rocs[[j]][i,2], rocs[[j]][i,3]), 1, 2)
        p<-rbind(p, pp)}
      }
      }
      if(nrow(p)>1)
      {
      med<-try(spatial.median(p[2:nrow(p),]))
      if(class(med)=="try-error")
        med<-c(median(p[2:nrow(p),1]), median(p[2:nrow(p),2]))
      data<-rbind(data, matrix(c(tempMax[i], med[1], med[2]), 1, 3))
      }
    }
    data
}

AUC<-function(ress){
  x<-ress[2:nrow(ress), 2]
  y<-ress[2:nrow(ress), 3]
  x2<-c(0.0, x[1:(length(x)-1)])
  y2<-c(1.0, y[1:(length(y)-1)])
  l<-(y+y2)*(x-x2)/2
  sum(l, na.rm=TRUE)
  }

PlotResults<-function(ress){
  N<-length(ress)
  for(i in 1:N)
  {
    x<-ress[[i]][2:nrow(ress[[i]]),2]
    y<-ress[[i]][2:nrow(ress[[i]]),3]
    ind<-which(!is.na(ROCBG[[1]])[,3] & !is.na(ROCBG[[1]])[,2])
    x<-x[ind]
    y<-y[ind]
    sortI<-sort(x, index.return=TRUE)
    plot(x, y, type="b", xlim=c(0, 1), ylim=c(0,1), xlab="sensitivity", ylab="specificity") ##, main=DecodeNoiseModels(2)TODO: add model type in main label
  }
}

Filzmoser<-function(num, bagplotAlpha=2.7, boxplotAlpha=2.7, bmethod = ""){
  epsSeq<-seq(from=0.05, to=0.5, length.out=10)
  data<-list()
  for(i in 1:10){
    #d<-GenerateContaminatedGauss(c(0,0), c(1, 1), 0.9, c(0,0), c(3, 3), -0.9, epsSeq[i], num) #10
    #d<-GenerateContaminatedGauss(c(0,0), c(1, 1), 0, c(0,0), c(3, 3), 0, epsSeq[i], num, fixedEps=TRUE) ##1
    d<-try(GenerateContaminatedGauss(c(0,0), c(1, 1), 0, c(3,0), c(1, 1), 0, epsSeq[i], num, fixedEps=TRUE)) ##1
    data[i]<-list(d)
    p1 <- proc.time()
    Sys.sleep(0.1)
    #print(proc.time() - p1)
    }
   # fix(data)
  res<-matrix(,ncol=5)
    dlength<- length(data)
  for(i in 1:dlength){
    df<-data[[i]]
    gc<-try(GetConfusionMatrix(df, boxplotAlpha))
    
    
    bgConf <- c()
    if(bmethod == "")
    {
      bp<-try(Bagplot(df, bagplotAlpha))
      bgConf<-try(GetOutliers(bp, df))##(TN, FN, FP, TP, sensitivity, specificity)
    }
    else
      bgConf <- try(GetConfusionMatrix(df, boxplotAlpha, method = bmethod))
   
    val1<- gc[6]
    val1W<-gc[2]/(gc[4]+gc[2])
    val2<- bgConf[6]
    val2W<-bgConf[2]/(bgConf[4]+bgConf[2])
    #print(c(val1, val1W, val2, val2W))
    res<-rbind(res, matrix(c(epsSeq[i], val1, val1W, val2, val2W), ncol=5))
  }
  res
}

IterateFilzmoser<-function(size)
{
  bpAlpha<-2.7
  #while maximum of wrong differences > 0.01
  i<-0
  eps<-1
  res<-NULL
  repeat{
    dmatrix<-Filzmoser(size, boxplotAlpha=bpAlpha)
    res1<-dmatrix[2:11,3]
    res2<-dmatrix[2:11,5]
    cond<-max(abs(res1-res2))
    if(cond<eps)
    {
      eps<-cond
      res<-list(table=dmatrix, alpha=bpAlpha)
    }
    if(i>10){
      break
      }
    bpAlpha<-bpAlpha+0.1
    i<-i+1
  }
  return(res)
}

PlotFilzmoser<-function(dmatrix){
  plot(dmatrix[1:10,1]*100, dmatrix[1:10,2]*100, type="l", ylim=c(0, 100), , xlab="Percent of contamination", ylab="Percent of outliers detected", col="red")
  lines(dmatrix[1:10,1]*100, dmatrix[1:10,4]*100, lty="dashed", col="blue") 
  lines(dmatrix[1:10,1]*100, dmatrix[1:10,3]*100, lty="dotdash", col="yellow")
  lines(dmatrix[1:10,1]*100, dmatrix[1:10,5]*100, lty="dotted", col="green")
}


ExploreFilzmoserPoints<-function()
{
  res <- matrix(rep(NA, 2*1000), ncol = 2)
  for(i in 1:1000)
  {
    d<-try(GenerateContaminatedGauss(c(0,0), c(1, 1), 0, c(3,0), c(1, 1), 0, 0.1, 50, fixedEps=TRUE)) ##1
    gc<-try(GetConfusionMatrix(d, 2.7))
    
    #bp<-try(Bagplot(d, 2.6))
    #gc<-try(GetOutliers(bp, d))
    if(class(GetROCPoints)!="try-error")
    {
      #if(length(gc) < 6)
       # print(gc)
    res[i, 1] <- gc[6]
    res[i, 2] <-gc[2]/(gc[4]+gc[2])
    }
    
  }
  return(res)
}

FilzmoserN <- function(SampleSize, NumPoints, bagplotAlpha=2.7, boxplotAlpha=2.3, bmethod = "", contType = "Center")
{
  axisPoints <- 1
  epsSeq<-seq(from=0.05, to=0.5, length.out=axisPoints)
  results <- matrix(rep(NA, 4 * axisPoints * NumPoints), ncol = 4*axisPoints)
  for(i in 1:axisPoints){
    for(j in 1:NumPoints)
    {
    #d<-GenerateContaminatedGauss(c(0,0), c(1, 1), 0.9, c(0,0), c(3, 3), -0.9, epsSeq[i], num) #10
      d <-c()
      
      if(contType == "Center")
    d<-GenerateContaminatedGauss(c(0,0), c(1, 1), 0, c(0,0), c(3, 3), 0, epsSeq[i], SampleSize) ##1
      else
    d<-try(GenerateContaminatedGauss(c(0,0), c(1, 1), 0, c(3,0), c(1, 1), 0, epsSeq[i], SampleSize)) ##1
      
    #data[i]<-list(d)
    #p1 <- proc.time()
    #Sys.sleep(0.1)
    #print(proc.time() - p1)
  
    # fix(data)
    res<-matrix(,ncol=5)
    
    df<-d
    gc<-try(GetConfusionMatrix(df, boxplotAlpha))
    
      oldSeed <- .Random.seed
    bgConf <- c()
    if(bmethod == "")
    {
      bp<-try(Bagplot(df, bagplotAlpha))
      bgConf<-try(GetOutliers(bp, df))##(TN, FN, FP, TP, sensitivity, specificity)
    }
    else
      bgConf <- try(GetConfusionMatrix(df, boxplotAlpha, method = bmethod))
      
      assign(".Random.seed", oldSeed, envir=globalenv())  
      
    val1<- gc[6]
    val1W<-gc[2]/(gc[4]+gc[2])
    val2<- bgConf[6]
    val2W<-bgConf[2]/(bgConf[4]+bgConf[2])
    #print(c(val1, val1W, val2, val2W))
    #res<-rbind(res, matrix(c(epsSeq[i], val1, val1W, val2, val2W), ncol=5))
    results[j, (i -1)*4 + 1 ] <- val1
    results[j, (i -1)*4 + 2 ] <- val1W
    results[j, (i -1)*4 + 3 ] <- val2
    results[j, (i -1)*4 + 4 ] <- val2W
    }
  }
  return(results)
}

GetConfusionMatrix2<-function(X, Y, O, a=2.7, method="elliptic"){
   
  ROC<-matrix(c(0,0,0), 1, 3)
  
  boxPlot <-c()
  if(method == "elliptic")
    boxPlot <- BuildBoxplot(X, Y, alpha=a)
  else
    if(method == "mve")
      boxPlot <- mvemethod(matrix(c(X, Y), ncol = 2 ))
  else
    if(method == "project")
      boxPlot <- promethod(data.frame(X = X, Y = X))
  else
    if(method == "mgv")
      boxPlot <- mgvmethod(matrix(c( X, Y), ncol = 2, byrow=TRUE) )
  else
    if(method == "mcd")
      boxPlot <- mcdmethod(data.frame(X = X, Y = Y))
  
  
  TN<- length( which( boxPlot$OutlierInds %in% which(O == 1) ) )
  FN<- length( which( boxPlot$OutlierInds %in% which(O == 0) ) )
  
  
  FP <- length(which(O==1)) - TN
  TP <- length(which(O==0)) - FN
  sensitivity<-NA
  if((TP+FN)>0)
    sensitivity<-TP/(TP+FN)
  specificity<-NA
  if((FP+TN)>0)
    specificity<-TN/(FP+TN)
  else
    if(FP+TN == 0 && FP >= 0 && TN >= 0)
      specificity <- 1
  #print(c(TN, FN, FP, TP))
  param<-c(TN, FN, FP, TP, sensitivity, specificity)
  return(param)
}

GetOutliers2<-function(bgplot, X, Y, O){
  bgOut<-try(bgplot$pxy.outlier)
  if(class(bgplot)!="bagplot")
  {
    fix(bgplot)
  }
  
  outInds <- which(O == 1)
  inInds <- which(O == 0)
  Out <- matrix(c(X[outInds], Y[outInds]), ncol=2)
  In <- matrix(c(X[inInds], Y[inInds]), ncol=2)
  
  TN<-0
  FN<-0
  
  if((length(bgOut)>0) && (nrow(Out)>0)){
    for(j in 1:nrow(bgOut)){
      for(i in 1:nrow(Out)){
        if(EqualPoints(Out[i,1], Out[i,2], bgOut[j,1], bgOut[j,2]))
          TN<-TN+1
      }
      
      for(i in 1:nrow(In)){
        if(EqualPoints(In[i,1], In[i,2], bgOut[j,1], bgOut[j,2]))
          FN<-FN+1
      }
    }
  }
  
  nbgOut<-0
  nOut<-0
  if(!is.null(nrow(bgOut)))
    nbgOut<-nrow(bgOut)
  
  if(!is.null(nrow(Out))) nOut<-nrow(Out)
  
  FP<-nOut - TN 
  TP<-nrow(In)-FN
  sensitivity<-NA
  if((TP+FN)>0)
    sensitivity<-TP/(TP+FN)
  specificity<-NA
  if((FP+TN)>0)
    specificity<-TN/(FP+TN)
  else
    if(FP+TN == 0 && FP >= 0 && TN >= 0)
      specificity <- 1
  param<-c(TN, FN, FP, TP, sensitivity, specificity)
  return(param)
}


FilzmoserNData<-function(data, bagplotAlpha=2.7, boxplotAlpha=2.3, bmethod = "ellipse")
{
  ddim <- dim(data)
  stopifnot( length(ddim) == 4)
  
  indxs <- matrix(1:ddim[1], ncol=1)
  
  tableX <- c()
  tableY <- c()
  tableO <- c()
  result <- matrix(, 0, 2)
  afun <- function(x){
    retVal <- c()
    xval <- tableX[x,]
    yval <- tableY[x,]
    oval <- tableO[x,]
    if(bmethod == "bagplot")
    {
      bp<-try(compute.bagplot(xval, yval, factor=bagplotAlpha, approx.limit=100000))
      retVal <- GetOutliers2(bp, xval, yval, oval)
    }
      else
      {
        retVal <- try(GetConfusionMatrix2(xval, yval, oval, boxplotAlpha, method = bmethod))
      }
    return(retVal)
  }
  
  for(i in 1:ddim[3])
  {
    tableX <- data[, ,i, 1]
    tableY <- data[, ,i, 2]
    tableO <- data[, ,i, 3]
    
    #val1<- gc[6]
    #val1W<-gc[2]/(gc[4]+gc[2])
    
    tresult <- apply(indxs, 1, afun)
    
    tresult2 <- matrix(c(tresult[6,], tresult[2,]/(tresult[4,] + tresult[2,])), ncol=2)
    hist(tresult[6,])
    tresult <- matrix(c(mean(tresult2[,1], na.rm=T), mean(tresult2[,2], na.rm = T) ), ncol=2)
    
    result<-rbind(result, tresult)
  }
  return(result)
}


FilzmoserNDataParallel<-function(data0, bagplotAlpha=2.7, boxplotAlpha=2.3, bmethod = "")
{
  ddim <- dim(data0)
  stopifnot( length(ddim) == 4)
  
  indxs <- matrix(1:ddim[1], ncol=1)
  
  tableX <- c()
  tableY <- c()
  tableO <- c()
  result <- matrix(, 0, 2)
  afun <- function(x, tXX, tYY, tOO){
    retVal <- c()
    xval <- tXX[x,]
    yval <- tYY[x,]
    oval <- tOO[x,]
    if(bmethod == "bagplot")
    {
      bp<-try(compute.bagplot(xval, yval, factor=bagplotAlpha, approx.limit=100000))
      retVal <- GetOutliers2(bp, xval, yval, oval)
    }
    else
    {
      retVal <- try(GetConfusionMatrix2(xval, yval, oval, boxplotAlpha, method = bmethod))
    }
    return(retVal)
  }
  
  mainp <- function(i)
  {
    tableX <- data0[, ,i, 1]
    tableY <- data0[, ,i, 2]
    tableO <- data0[, ,i, 3]
    
    tresult <- apply(matrix(1:dim(tableX)[1], ncol=1), 1, afun, tXX=tableX, tYY=tableY, tOO=tableO)
    
    
    tresult2 <- matrix(c(tresult[6,], tresult[2,]/(tresult[4,] + tresult[2,])), ncol=2)
   
    tresult <- matrix(c(mean(tresult2[,1], na.rm=T), mean(tresult2[,2], na.rm = T) ), ncol=2)
    if(dim(tresult)[2] != 2)
      tresult <- matrix(numeric(0), ncol=2)
    #result<-rbind(result, tresult)
    return (tresult)
  }
  
  cores <- 3
  
  cl <- makeCluster(cores)
  assign("data0", data0, pos=1)
  assign("tableX", tableX, pos=1)
  assign("tableY", tableY, pos=1)
  assign("tableO", tableO, pos=1)
  assign("boundOpt", boundOpt, pos=1)
  assign("optimize", optimize, pos=1)
  assign("indxs", indxs, pos=1)
  assign("GetConfusionMatrix2", GetConfusionMatrix2, pos=1)
  assign("BuildBoxplot", BuildBoxplot, pos=1)
  assign("EstimateLocation", EstimateLocation, pos=1)
  assign("CorrFQn", CorrFQn, pos=1)
  assign("M_1_SQRT2", M_1_SQRT2, pos=1)
  assign("FQn", FQn, pos=1)
  assign("EstimateScale", EstimateScale, pos=1)
  assign("EstimateDistance", EstimateDistance, pos=1)
  assign("EstimateCovMatr", EstimateCovMatr, pos=1)
  assign("GetHalfPointsFromCenterLocation", GetHalfPointsFromCenterLocation, pos=1)
  assign("GetFencePoints", GetFencePoints, pos=1)
  assign("mgvmethod", mgvmethod, pos=1)
  assign("outmgv", outmgv, pos=1)
  assign("mgvar", mgvar, pos=1)
  assign("elimna", elimna, pos=1)
  assign("apgdis", apgdis, pos=1)
  assign("mean", mean, pos=1)
  assign("median", median, pos=1)
  assign("mad", mad, pos=1)
  assign("outer", outer, pos=1)
  assign("outbox", outbox, pos=1)
  assign("order", order, pos = 1)
  assign("qchisq", qchisq, pos=1)
  assign("sqrt", sqrt, pos = 1)
  assign("rep", rep, pos=1)
  clusterExport(cl, c("data0", "tableX", "tableY", "tableO", "GetConfusionMatrix2", "BuildBoxplot", "CorrFQn", "EstimateLocation", 
                      "FQn", "M_1_SQRT2", "EstimateScale", "EstimateDistance", "EstimateCovMatr", "GetHalfPointsFromCenterLocation", 
                      "GetFencePoints", "mgvmethod", "outmgv", "mgvar", "elimna", "apgdis", "mean", "median", "mad", "outer", "outbox", "order",
                      "qchisq", "sqrt", "rep"))
  clusterEvalQ(cl, require("stats"))
  clusterEvalQ(cl, source("Rallfun-v20.R"))
  result<-clusterApplyLB(cl, seq(1, ddim[3], by=2), mainp)
  stopCluster(cl)
  
  return(result)
}


FindOptimalBounds<-function(data0)
{
  ddim <- dim(data0)
  stopifnot( length(ddim) == 4)
  
  indxs <- matrix(1:ddim[1], ncol=1)
  bmethod <- "elliptic"
  
  tableX <- c()
  tableY <- c()
  tableO <- c()
  result <- 0
  cores <- 2
  
  boundOpt <- function(kBound, tX=tableX, tY=tableY, tO=tableO)
  {
    afun <- function(x, tXX, tYY, tOO){
      retVal <- c()
      xval <- tXX[x,]
      yval <- tYY[x,]
      oval <- tOO[x,]
      retVal <- try(GetConfusionMatrix2(xval, yval, oval, kBound, method = bmethod))
      return(retVal)
    }
    
    sinds <- sample(nrow(tX), 100)
    
    tresult <- apply(matrix(sinds, ncol=1), 1, afun, tXX=tX, tYY=tY, tOO=tO)
    #return(tresult)
    tresult2 <- matrix(c(tresult[6,], tresult[5,]), ncol=2)#tresult[2,]/(tresult[4,] + tresult[2,])
    sens<- mean(tresult2[,1], na.rm=T)
    spec<- mean(tresult2[,2], na.rm=T)
    
    HResult <- 2 * sens * spec / (sens + spec)
    return(HResult)
  }
  
  #for(i in 1:ddim[3])
    mainp <- function(i)
  {
    
    tableX <- data0[, ,i, 1]
    tableY <- data0[, ,i, 2]
    tableO <- data0[, ,i, 3]
    
    #val1<- gc[6]
    #val1W<-gc[2]/(gc[4]+gc[2])
    resultH <- optimize(f=boundOpt, interval = c(0.5, 4), tX=tableX, tY=tableY, tO=tableO,  maximum=T)
    #resultH <- boundOpt(1, tX=tableX, tY=tableY, tO=tableO)
    #print(resultH)
    #result<-rbind(result, matrix( c(resultH$objective, resultH$maximum), ncol = 2) )
    #return (resultH)
    return(resultH)
  }
  
  
  
  cl <- makeCluster(cores)
  assign("data0", data0, pos=1)
  assign("tableX", tableX, pos=1)
  assign("tableY", tableY, pos=1)
  assign("tableO", tableO, pos=1)
  assign("boundOpt", boundOpt, pos=1)
  assign("optimize", optimize, pos=1)
  assign("indxs", indxs, pos=1)
  assign("GetConfusionMatrix2", GetConfusionMatrix2, pos=1)
  assign("BuildBoxplot", BuildBoxplot, pos=1)
  assign("EstimateLocation", EstimateLocation, pos=1)
  assign("CorrFQn", CorrFQn, pos=1)
  assign("M_1_SQRT2", M_1_SQRT2, pos=1)
  assign("FQn", FQn, pos=1)
  assign("EstimateScale", EstimateScale, pos=1)
  assign("EstimateDistance", EstimateDistance, pos=1)
  assign("EstimateCovMatr", EstimateCovMatr, pos=1)
  assign("GetHalfPointsFromCenterLocation", GetHalfPointsFromCenterLocation, pos=1)
  assign("GetFencePoints", GetFencePoints, pos=1)
  clusterExport(cl, c("data0", "tableX", "tableY", "tableO", "GetConfusionMatrix2", "BuildBoxplot", "CorrFQn", "EstimateLocation", "FQn", "M_1_SQRT2", "EstimateScale", "EstimateDistance", "EstimateCovMatr", "GetHalfPointsFromCenterLocation", "GetFencePoints"))
  #clusterEvalQ(cl, {library(fastqn, pcaPP, ICSNP)})
  result<-clusterApply(cl, 1:ddim[3], mainp)
  stopCluster(cl)
  
  return(result)
}


TranslateFN <- function(fz)
{
  t <- apply(fz, 2, mean)
  tm <- matrix(t, ncol = 4, byrow = TRUE)
  #PlotFilzmoser(tm)
  tm <- cbind(matrix(seq(from = 0.05, to = 0.5, length.out = 10), ncol=1), tm)
  tm <- rbind(matrix(rep(NA, 5), ncol=5), tm)
  
  return(tm)
}

ComputeFMeasure<-function(fz)
{
  retVal <- fz[,1:3]
  retVal[,2]<- 2 * fz[,2] * (1-fz[,3])/(fz[,2] - fz[,3] + 1)
  retVal[,3]<- 2 * fz[,4] * (1-fz[,5])/(fz[,4] - fz[,5] + 1)
  return(retVal)
}

plotBounds<-function()
{
  library("ggplot2")
  m100 <- matrix(unlist(bounds.DataS100.N10000), ncol=2, byrow=T)
  m100Shift <- matrix(unlist(bounds.DataS100.N10000Shift), ncol=2, byrow=T)
  m50 <- matrix(unlist(bounds.DataS50.N10000), ncol=2, byrow=T)
  m50Shift <- matrix(unlist(bounds.DataS50.N10000Shift), ncol=2, byrow=T)
  m1000 <- matrix(unlist(bounds.DataS1000.N10000_1), ncol=2, byrow=T)
  m1000Shift <- matrix(unlist(bounds.DataS1000.N10000_1Shift), ncol=2, byrow=T)
  
  dtFrame <- data.frame(
    H = c(m50[,2], m50Shift[,2], m100[,2], m100Shift[,2], m1000[,2], m1000Shift[,2]),
    T0 = c(m50[,1], m50Shift[,1], m100[,1], m100Shift[,1], m1000[,1], m1000Shift[,1]),
    Contamination = rep(seq(from=0.05, to=0.5, length.out=10), times=6),
    CType = rep(rep(c("Scale", "Shift"), each = 10 ), times = 3, length.out=60),
    Size = rep(c(50, 100, 1000), each=20, length.out=60)
    ) 
  return(dtFrame)
  #ggplot(data = dtFrame) + geom_bar(aes(x=Contamination, y = H, fill=CType), stat="identity", position=position_dodge()) +theme(axis.title.x = element_text(face="bold", colour="#990000", size=26), axis.text.x  = element_text(size=20), axis.title.y = element_text(face="bold", colour="#990000", size=26), axis.text.y  = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))+ geom_line(aes(x=Contamination, y = T0, color=CType)) + facet_grid(.~Size)
}

combineParallel <- function(ctype = "center")
{
  par1 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_1), ncol=2, byrow=T)
  par2 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_2), ncol=2, byrow=T)
  par3 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_3), ncol=2, byrow=T)
  par4 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_4), ncol=2, byrow=T)
  par5 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_5), ncol=2, byrow=T)
  par6 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_6), ncol=2, byrow=T)
  par7 <- matrix(unlist(resultparallel_DataS1000.N10000Shift_7), ncol=2, byrow=T)
  
  par1C <- matrix(unlist(resultparallel_DataS1000.N10000_1), ncol=2, byrow=T)
  par2C <- matrix(unlist(resultparallel_DataS1000.N10000_2), ncol=2, byrow=T)
  par3C <- matrix(unlist(resultparallel_DataS1000.N10000_3), ncol=2, byrow=T)
  par4C <- matrix(unlist(resultparallel_DataS1000.N10000_4), ncol=2, byrow=T)
  par5C <- matrix(unlist(resultparallel_DataS1000.N10000_5), ncol=2, byrow=T)
  par6C <- matrix(unlist(resultparallel_DataS1000.N10000_7), ncol=2, byrow=T)
  
  parMean1 <- rowMeans(cbind(par1[,1], par2[,1], par3[,1], par4[,1], par5[,1], par6[,1], par7[,1]))
  parMean2 <- rowMeans(cbind(par1[,2], par2[,2], par3[,2], par4[,2], par5[,2], par6[,2], par7[,2]))
  
  parMean1C <- rowMeans(cbind(par1C[,1], par2C[,1], par3C[,1], par4C[,1], par5C[,1], par6C[,1]))
  parMean2C <- rowMeans(cbind(par1C[,2], par2C[,2], par3C[,2], par4C[,2], par5C[,2], par6C[,2]))
  
  if(ctype == "center")
    return(2 * parMean1C *(1 - parMean2C) / (parMean1C + 1 - parMean2C))
  
  return(2 * parMean1 * (1 - parMean2) / (parMean1 + 1 - parMean2))
}

plotComparison <- function(c1_50, c1_100, c1_1000, c2_50, c2_100, c2_1000, c1_name, c2_name )
{
  mc1_50 <- c1_50
  mc1_100 <- c1_100
  mc1_1000 <- c1_1000
  mc2_50 <- c2_50
  mc2_100 <- c2_100
  mc2_1000 <- c2_1000
  
  if(is.matrix(c1_50) && dim(c1_50)[2] == 2)
    mc1_50 <- 2 * c1_50[,1] * (1- c1_50[,2]) / (c1_50[,1] + 1 - c1_50[,2])    
  if(is.matrix(c1_100) && dim(c1_100)[2] == 2)
    mc1_100 <- 2 * c1_100[,1] * (1- c1_100[,2]) / (c1_100[,1] + 1 - c1_100[,2])
  if(is.matrix(c1_1000) && dim(c1_1000)[2] == 2)
    mc1_1000 <- 2 * c1_1000[,1] * (1- c1_1000[,2]) / (c1_1000[,1] + 1 - c1_1000[,2])
  
  if(is.matrix(c2_50) && dim(c2_50)[2] == 2)
    mc2_50 <- 2 * c2_50[,1] * (1- c2_50[,2]) / (c2_50[,1] + 1 - c2_50[,2])    
  if(is.matrix(c2_100) && dim(c2_100)[2] == 2)
    mc2_100 <- 2 * c2_100[,1] * (1- c2_100[,2]) / (c2_100[,1] + 1 - c2_100[,2])
  if(is.matrix(c2_1000) && dim(c2_1000)[2] == 2)
    mc2_1000 <- 2 * c2_1000[,1] * (1- c2_1000[,2]) / (c2_1000[,1] + 1 - c2_1000[,2])
  
  dtFrame <- data.frame(
    H = c(mc1_50, mc1_100, mc1_1000, mc2_50, mc2_100, mc2_1000),
    Contamination = rep(seq(from=0.05, to=0.5, length.out=10), times=6),
    Boxplot = rep(c(c1_name, c2_name), each = 30),
    Size = rep(c(50, 100, 1000), each = 10, times = 2)
    )
  
  p <- ggplot(data = dtFrame) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=26), axis.text.x  = element_text(size=20), axis.title.y = element_text(face="bold", colour="#990000", size=26), axis.text.y  = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))+ geom_line(aes(x=Contamination, y = H, color=Boxplot)) + facet_grid(.~Size)
}