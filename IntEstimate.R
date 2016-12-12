##pivec<-seq(from=0, to=6.282, length.out=1000)
REllipse<-function(center, correl, sigma1, sigma2, piseq, Phi=NA){
  phi<-0.0
  if(!is.na(Phi)){
    phi<-Phi
    }else
  if(correl != 0.0){
    aphi<-2*correl*(sigma1*sigma2)/(sigma1^2-sigma2^2)
    phi<-atan(aphi)/2
  }
  
  A<-sigma1
  B<-sigma2
  RE<-A*B/sqrt((B*cos(piseq-phi))^2+(A*sin(piseq- phi))^2)
  
  return(RE)
}

Epoints<-function(sigma1, sigma2, piseq)
{
  phi<-0.2449787
  x<-sigma1*cos(phi)*cos(piseq-phi)-sigma2*sin(phi)*sin(piseq-phi)
  y<-sigma1*sin(phi)*cos(piseq-phi)+sigma2*cos(phi)*sin(piseq-phi)
  retVal<-matrix(c(x,y), ncol=2)
  return(retVal)
  }

FindIntersect<-function(a, y0, x0, tgPhi){
  if(is.infinite(a)){
      val<-c(x0, x0*tgPhi)
      return(val)
    }
    else
    {
      A<-matrix(c(1, 1, a, -tgPhi), ncol=2)
      b<-c(y0+a*x0, 0)
      res<-try(solve(A, b))
      if(class(res)=="try-error")
        print(c(a, y0, x0, tgPhi))
      val<-c(res[2], res[1])
      return(val)
    }
}

RHull<-function(center, hullPts, piseq){
  #print(hullPts)
  hullPts<-matrix(c(hullPts[,1] - center[1], hullPts[,2] - center[2]), ncol=2)
  hIndex<-c(nrow(hullPts), 1:(nrow(hullPts)-1))
  
  k<- -(hullPts[,2]-hullPts[hIndex,2])/(hullPts[,1]-hullPts[hIndex,1])
  #fix(k)
  tgPhi<- tan(piseq)
  PhiRange<- (atan2(hullPts[,2], hullPts[,1]) + 2*pi)%%(2*pi)

  PhiSorted<- sort(PhiRange, index.return=TRUE)
  #fix(PhiSorted)
  i<-1
  #print(matrix(c(hIndex, hullPts[hIndex,1], hullPts[hIndex,2], k, PhiSorted$x), ncol=5))
  
  j<- 1##PhiSorted$ix[1]
  Phi<-c(PhiSorted$x,2*pi)
  M<-length(Phi)
  pts<-matrix(c(0,0), ncol=2)
  for(i in 1:length(piseq)){
    ##print(c(Phi[j], piseq[i]))
    while(Phi[j]<piseq[i]){
      j<-(j+1)%%(M + 1)
      if(j == 0) j<-1
    }
    
    ##print(c(k[j],hullPts[j, 2], hullPts[j, 1], tgPhi[i]))
    J<-if(j%%M== 0) 1 else j%%(M)
    
    IPt<-FindIntersect(k[PhiSorted$ix[J]],hullPts[PhiSorted$ix[J], 2], hullPts[PhiSorted$ix[J], 1], tgPhi[i])
    
    pts<-rbind(pts, matrix(c(IPt[1], IPt[2]), ncol=2))
    
  }
  
  sqrt(pts[2:nrow(pts),1]^2+pts[2:nrow(pts),2]^2)
  
}

IntegrateR<-function(RE, RH){
  if(length(RE)!= length(RH)) print("Different lengths")
  dif<-RH-RE
  #print(matrix(c(RH,RE, dif), ncol=3))
  I<-(dif/RE)^2
  dPhi<-2*pi/length(RE)
  I<-I*dPhi/2
  I<-sum(I, na.rm=TRUE)
  I/(2*pi)
}

BuildCircleHull<-function(Radius, dPhi){
  pivec<-seq(from=0.2, to=6.282, length.out=10)
  matrix(c(cos(pivec), sin(pivec)), ncol=2)
}

BuildEllipseHull<-function(R1, R2){
  piseq<-seq(from=0, to=6.24, length.out=100)
  phi<-pi/4
  theta<-piseq
  x<-R1*cos(phi)*cos(theta)-R2*sin(phi)*sin(theta)
  y<-R1*sin(phi)*cos(theta)+R2*cos(phi)*sin(theta)
  matrix(c(x, y), ncol=2)
}

TestIntegral<-function(){
  pivec<-seq(from=0, to=6.24, length.out=100)
  #CH<-BuildCircleHull(1, 0)
  
  CH<-BuildEllipseHull(2.0, 2.0)
  plot(CH, type="l")
  REt<-REllipse(c(0,0), 1.0, 2.0, 1.0, pivec, Phi=pi/4)
  #print(REt)
  RHt<-RHull(c(0, 0), CH, pivec)
  
  #lines(Epoints(2.0, 1/2, pivec))
  IntegrateR(REt, RHt)
}

CalculateIntegral<-function(EllipseParams, HingeHull){
  pivec<-seq(from=0, to=6.282, length.out=100)
  REt<-REllipse(EllipseParams$center, EllipseParams$corr, EllipseParams$A, EllipseParams$B, pivec)
  if(any(is.na(REt)))
     print(c(EllipseParams$center, EllipseParams$corr, EllipseParams$A, EllipseParams$B))
  RHt<-RHull(HingeHull$center, HingeHull$Pts, pivec)
  ##if(any(is.na(RHt)))
  ##  print(HingeHull)
  if(length(REt) != length(RHt))
    print(matrix(c(REt, RHt), ncol=2))
  IntegrateR(REt, RHt)
  }
  
RunMonteCarlo<-function(num, bptype="Boxplot"){
  
  fres<-matrix(c(0,0,0,0), ncol=4)
  #for all models of noise
   
  for(j in 1:11){
   
    print(j)
    meanV<-0
    meanVec<-c()
    res<-0
    resmed<-c()
    inputlist<-list()
    # add num lists of generated noise
    for(i in 1:num){
      SignalNoise<-GenerateNoiseByModel(type=j, N=500)
    
      x<-SignalNoise$X
	    y<-SignalNoise$Y
      inputlist[i]<-list(matrix(c(x,y), ncol=2))
    }
  # calculate boxplot for each list of data and retrieve results
    for(i in 1:num)
  {
    SigParams<-GetNoiseParametersByModel(type=j)
    ##list(CovariationMCD = CovMat, Phi=phi, Points=dataMatrix, Center=boxplotCenter, BoxplotScale=boxplotScale, Hinge=hinge, Fence = fence, Outlier=outliers, OutlierInds=outliersInds)
    bp<-list()
    x<-inputlist[[i]][,1]
    y<-inputlist[[i]][,2]
    ##print(inputlist[[i]])
    k<-2.7
    if(bptype=="Boxplot"){
      bp<-BuildBoxplot(x, y, k)
    }
    else
    {
      oldSeed <- .Random.seed 
      bp<-try(BG(x,y, al=k))
      # Push random seed to old state 
      assign(".Random.seed", oldSeed, envir=globalenv()) 
    }
    ##try(PlotBoxplot(bp, scatterPlot=FALSE))
    EP<-list(center=c(0,0), corr=SigParams$rro1, A=k*SigParams$sigma1[1], B=k*SigParams$sigma1[2])
    ##print(EP)
    ##print(matrix(c(x,y), ncol=2))
    IVal<-CalculateIntegral(EP, list(Pts=bp$Fence[(nrow(bp$Fence)-1):1,], center=bp$Center))
    ##print(IVal)
  dist<-sqrt(bp$Center[1]^2+bp$Center[2]^2)
  meanV<- meanV + dist
  if(is.na(meanV))
    print(c(dist, bp$center[1], bp$center[1]))
    meanVec<-c(meanVec, dist)
    res<-res+IVal
    resmed<-c(resmed, IVal)
  }
  ##print(c(res/num, median(resmed)))
  fres<-rbind(fres, matrix(c(res/num, median(resmed), meanV/num, median(meanVec)), ncol=4))
  }
  
  return(fres)
}

MonteCarloEllipseEval2 <- function(byPC=F)
{
  resultList<-vector(mode="list", length=11)
  for(j in 1:11)
  {
    NMod <- 100
    inputlist <- lapply(1:NMod, function(x){
      SignalNoise <- GenerateNoiseByModel(type=j, N=500)
      x<-SignalNoise$X
      y<-SignalNoise$Y
      return(matrix(c(x,y), ncol=2))
    })
    
    
    resultMat <- matrix(rep(0, NMod*4), ncol=4)
    
    i<-0
    for(dataM in inputlist)
    {
      i <- i + 1
      
      bp<-BuildBoxplot(dataM[,1], dataM[,2], 3)
      bag <- BG(dataM[,1], dataM[,2])
      rel <- relplot(dataM[,1], dataM[,2], plotit=F)
      quel <- quelplot(dataM[,1], dataM[,2], plotit=F)
      
      SigParams<-GetNoiseParametersByModel(type=j)
      
      EP<-list(center=c(0,0), corr=SigParams$rro1, A=3*SigParams$sigma1[1], B=3*SigParams$sigma1[2])
      
      bpVal<-CalculateIntegral(EP, list(Pts=bp$Fence[(nrow(bp$Fence)-1):1,], center=bp$Center))
      if(byPC)
      {
        #REllipse<-function(center, correl, sigma1, sigma2, piseq, Phi=NA){
        bpEl <- REllipse(c(0,0), bp$Correlation$Corr, 3*bp$BoxplotScale[1], 3*bp$BoxplotScale[2], seq(from=0, to=6.282, length.out=100))
        REt<-REllipse(c(0,0), SigParams$rro1, 3*SigParams$sigma1[1], 3*SigParams$sigma1[2], seq(from=0, to=6.282, length.out=100))
        bpVal<- IntegrateR(REt, bpEl)
      }
      
      bagVal<-CalculateIntegral(EP, list(Pts=bag$Fence[(nrow(bag$Fence)-1):1,], center=bag$Center))
      
      relVal <-CalculateIntegral(EP, list(Pts=rel$fence[(nrow(rel$fence)-1):1,], center=rel$Center))
      quelVal <-CalculateIntegral(EP, list(Pts=quel$fence[(nrow(quel$fence)-1):1,], center=quel$Center))
      resultMat[i, ] <- c(bpVal, bagVal, relVal, quelVal)
      
    }
    print(mean(resultMat[,1]))
    resultList[j] <- list(tab = resultMat)
    
  }
  return(resultList)
}

BG<-function(x, y, al=3)
{
  oldSeed <- .Random.seed
  bg<-compute.bagplot(x, y, factor=al)
  # Push random seed to old state 
  assign(".Random.seed", oldSeed, envir=globalenv()) 
  
  list(Center=bg$center, Hinge=rbind(bg$hull.bag, matrix(c(0,0),ncol=2)), Fence=rbind(bg$hull.loop, matrix(c(0,0),ncol=2)))
}
