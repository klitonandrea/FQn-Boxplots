promethod<-function(x)
{
  res <- outpro(x, plotit = FALSE)
  return (list(Points = x, OutlierInds=res$out.id))
}

mgvmethod<-function(x)
{
  res <- outmgv(x, plotit = FALSE)
  return(list(Points = x, OutlierInds = res$out.id))
}

mvemethod <- function(x)
{
  res<-outmve(x, plotit = FALSE)
  return(list(Points = x, OutlierInds = res$out.id))
}

mcdmethod <- function(x)
{
  res <- out(x, cov.fun=cov.mcd, plotit = FALSE)
  return(list(Points = x, OutlierInds = res$out.id))
}