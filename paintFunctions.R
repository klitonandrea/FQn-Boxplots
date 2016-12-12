paintGaussian <- function()
{
  dat <- dnorm(1000)
  curve(dnorm, from = -4, to = 4, n = 1000, axes = T, xlab = "", ylab = "")
  dat25 <- qnorm(0.25)
  dat75 <- qnorm(0.75)
  
  x <- seq(from=dat25, to=dat75, length = 500)
  polygon(c(dat25, x, dat75), c(0, dnorm(x), 0))#, col ="red")
  
  lx <- seq(from = -4, to = dat25, length = 500)
  polygon(c(-4, lx, dat25), c(0, dnorm(lx), 0))#, col = "blue")
  
  ux <- seq(from = dat75, to = 4, length = 500)
  polygon(c(dat75, ux, 4), c(0, dnorm(ux), 0))#, col = "blue")
  boxplot(rnorm(1000), add=TRUE, horizontal = TRUE, axes = F)
}

paintHM <- function()
{
  library("rgl")
  H <- function(r, p){ ifelse(r ==0 && p==0, 0, 2 * r * p / (r + p))}
  recall <- seq(from = 0, to= 1, length=10)
  precision <- seq(from = 0, to = 1, length = 10)
  f <- matrix(rep(0, length(recall)*length(precision)), ncol = length(recall))
  
  for( i in 1:length(recall))
    for(j in 1:length(precision))
      f[j, i] <- H(recall[i], precision[j])
  
  zlim <- range(precision * 100)
  zlen <- zlim[2] - zlim[1] + 1
  colorlut <- terrain.colors(zlen)
  col <- colorlut[ f * 100 - zlim[1] + 1]
  persp(recall, precision, f, zlab = "H", xlab = "Specificity", ylab = "Sensitivity", phi = 30, theta = -30, ticktype = "detailed", col = col, cex.lab = 1.5, cex.axis=1.)
  
  ##s <- surfaceTriangles(recall, precision, f)
  ##drawScene(s)
  
  #library("rgl")
  #open3d()
  #surface3d(recall, precision, f, color = col, back = "lines")
}


ScaleCResults<-data.frame(
  N = rep(factor(c(20, 50, 100, 1000, 10000)), each = 4),
  Boxplot = rep(c("Tukey", "MAD", "FQn", "Grubbs"), 5),
  H = c(0.64, 0.67, 0.66, 0.17, 0.72, 0.72, 0.72, 0.29, 0.72, 0.73, 0.72, 0.3, 0.72, 0.73, 0.72, 0.3, 0.72, 0.73, 0.73, 0.3)
)

ShiftCResults<-data.frame(
  N = rep(factor(c(20, 50, 100, 1000, 10000)), each = 4),
  Boxplot = rep(c("Tukey", "MAD", "FQn", "Grubbs"), 5),
  H = c(0.75, 0.73, 0.73, 0.32, 0.79, 0.8, 0.79, 0.39, 0.8, 0.8, 0.81, 0.4, 0.8, 0.8, 0.81, 0.39, 0.8, 0.8, 0.81, 0.39)
)

ShiftEpsResults <- data.frame(
  eps = rep(c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), each = 4),
  Boxplot = rep(c("Tukey", "MAD", "FQn", "Grubbs"), 6),
  H = c(0.63, 0.65, 0.67, 0.65, 0.62, 0.65, 0.67, 0.56, 0.59, 0.6, 0.61, 0.41, 0.55, 0.56, 0.56, 0.31, 0.51, 0.52, 0.5, 0.25, 0.43, 0.44, 0.4, 0.21)
)

RealResults<-data.frame(
  Measurement = rep(c("%idle", "tpc","bread/s+bwrtn/s", "rxpck/s+txpck/s"), each=3),
  Boxplot = rep(c("Tukey", "MAD", "FQn"), 4),
  H = c(0.51, 0.55, 0.58, 0.47, 0.57, 0.57, 0.47, 0.56, 0.56, 0.2, 0.33, 0.31)
)


ResultPaint <- function(Rdata)
{
  library("ggplot2")
  p <- ggplot(data=Rdata, aes(x=N, y=H, fill=Boxplot)) + geom_bar(stat="identity", position=position_dodge())
  return(p)
}

PaintContaminated<-function()
{
  cont<-function(x)(dnorm((x-3)/3 )/3)
  curve(dnorm, from = -9, to = 9, n = 1000, axes = T, xlab = "", ylab = "", col='red', lwd=3)
  curve(cont, from = -9, to = 9, n = 1000, axes = T, xlab = "", ylab = "", col='blue', lwd=3, add=T)
}