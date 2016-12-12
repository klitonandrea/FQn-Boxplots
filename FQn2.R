##############################################################################
EllipsePrincAxis <- function(x, y, r, m1, m2, sig1, sig2){
  tan2phi <- 2 * r * sig1 * sig2 / ((sig1^2) - (sig2^2))
  phi <- atan(tan2phi)
  res <- (x - m1)*cos(phi) + (y - m2)*sin(phi)
  res
}
##############################################################################
EllipseConjAxis <- function(x, y, r, m1, m2, sig1, sig2){
  tan2phi <- 2 * r * sig1 * sig2 / ((sig1^2) - (sig2^2))
  phi <- atan(tan2phi)
  res <- -(x - m1)*sin(phi) + (y - m2)*cos(phi)
  res
}
##############################################################################
MakeStack<-function(x){
  pop<-function() if (length(x)>1) x<<-x[-1]
  push<-function(y) x<<-c(y,x)
  value<-function() x
  rval<-list(pop=pop,push=push,value=value)
  class(rval)<-"myStack"
  rval
}
##############################################################################
IsLeft <- function(p0, p1, p2){
  is_left <- (p1[1] - p0[1])*(p2[2] - p0[2]) - (p2[1] - p0[1])*(p1[2] - p0[2])
  if(is_left >= 0){
    FALSE
  }
  else{
    TRUE
  }
}
#########################################################################################
GetPolarAngle <- function(p0, p){
  if(identical(p0, p)){
    a <- 0
  }
  else{
    tg_a <- (p[2] - p0[2]) / (p[1] - p0[1])
    a <- atan(tg_a)
    a <- a * 180 / pi
    if(p[1] < p0[1]){
      a <- 180 + a
    }
  }
  a
}
##############################################################################
Graham <- function(Data) {
  #bottom point (p0)
  Data_sort <- Data[order(Data$v2, Data$v1), ]
  imax <- length(Data_sort$v1)
  38
  p0 <- c(Data_sort$v1[1], Data_sort$v2[1])
  Angles <- vector(mode = "numeric", length = imax)
  for (i in 1:imax){
    p <- c(Data$v1[i], Data$v2[i])
    Angles[i] <- GetPolarAngle(p0, p)
  }
  #angle sort
  Data_struct <- data.frame(v1 = Data$v1,
                            v2 = Data$v2,
                            v3 = Angles)
  Data_struct <- Data_struct[order(Data_struct$v3),]
  #stack contains only indexes
  stack <- MakeStack(imax)
  stack$push(1)
  stack$value()
  i <- 2
  while(i <= imax){
    l <- length(stack$value())
    j <- stack$value()[1]
    k <- stack$value()[2]
    p_cur <- c(Data_struct$v1[i], Data_struct$v2[i])
    p1 <- c(Data_struct$v1[j], Data_struct$v2[j])
    p2 <- c(Data_struct$v1[k], Data_struct$v2[k])
    if(IsLeft(p_cur, p1, p2)){
      stack$push(i)
      i <- i + 1
    }
    else{
      stack$pop()
    }
  }
  res_Data <- data.frame(v1 = Data_struct$v1[stack$value()],
                         v2 = Data_struct$v2[stack$value()])
  res_Data
}
##############################################################################
Qn_boxplot <- function(Data, draw_points)
{
  Qn_x <- FQn(Data$v1)
  Qn_y <- FQn(Data$v2)
  #the final robust estimate of the covariance matrix
  # and correlation matrix
  rob_estim <- covPairs(Data)
  cov <- matrix(c(rob_estim$cov), ncol = 2, byrow = TRUE)
  cor <- matrix(c(rob_estim$cor), ncol = 2, byrow = TRUE)
  center <- rob_estim$center
  dist <- sqrt(rob_estim$dist)
  #inner
  imax <- length(Data$v1)
  D_dist <- data.frame(v1 = dist,
                       v2 = Data$v1,
                       v3 = Data$v2)
  D_dist <- D_dist[order(D_dist$v1, D_dist$v2, D_dist$v3), ]
  jmax <- imax %/% 2
  inner <- data.frame(v1 = D_dist$v2[1:jmax],
                      v2 = D_dist$v3[1:jmax])
  u <- EllipsePrincAxis(Data$v1, Data$v2, cor[2], center[1], center[2], Qn_x, Qn_y)
  v <- EllipseConjAxis(Data$v1, Data$v2, cor[2], center[1], center[2], Qn_x, Qn_y)
  med_u <- median(u)
  med_v <- median(v)
  u <- sort(u)
  v <- sort(v)
  Qn_u <- FQn(u)
  Qn_v <- FQn(v)
  x_L <- max(u[1], med_u - 3 * Qn_u)
  y_L <- max(v[1], med_v - 3 * Qn_v)
  x_R <- min(u[imax], med_u + 3 * Qn_u)
  y_R <- min(v[imax], med_v + 3 * Qn_v)
  #outer
  outer_indexes <- vector(mode = "numeric", length = 0)
  outlier_indexes <- vector(mode = "numeric", length = 0)
  jmax <- jmax + 1
  for (i in jmax : imax){
    p_x <- EllipsePrincAxis(D_dist$v2[i], D_dist$v3[i], cor[2],
                            center[1], center[2], Qn_x, Qn_y)
    p_y <- EllipseConjAxis(D_dist$v2[i], D_dist$v3[i], cor[2],
                           center[1], center[2], Qn_x, Qn_y)
    if((p_x < x_R)&&(p_x > x_L)&&(p_y < y_R)&&(p_y > y_L))
      outer_indexes <- c(outer_indexes, i)
    else
      outlier_indexes <- c(outlier_indexes, i)
  }
  outer <- data.frame(v1 = D_dist$v2[outer_indexes],
                      v2 = D_dist$v3[outer_indexes])
  outliers <- data.frame(v1 = D_dist$v2[outlier_indexes],
                         v2 = D_dist$v3[outlier_indexes])
  #convex hull(Graham's scan)
  inner_con_hull <- Graham(inner)
  in_outer <- data.frame(v1 = c(inner$v1, outer$v1),
                         v2 = c(inner$v2, outer$v2))
  outer_con_hull <- Graham(in_outer)
  #plot
  plot(NULL, xlab="X", ylab="Y",
       xlim = c(min(outer_con_hull$v1, outliers$v1), max(outer_con_hull$v1, outliers$v1)),
       ylim = c(min(outer_con_hull$v2, outliers$v2), max(outer_con_hull$v2, outliers$v2)))
  polygon(outer_con_hull$v1, outer_con_hull$v2, border = 'gray82', col = 'gray82')
  polygon(inner_con_hull$v1, inner_con_hull$v2, border = 'gray57', col = 'gray57')
  points(outliers, pch = 8, col = 2)
  points(center[1], center[2], pch = 3, col = 1)
  if(draw_points == TRUE)
  {
    points(outer, pch = 20)
    points(inner, pch = 20)
  }
}
