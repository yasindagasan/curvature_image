# R function
###############################################################################
#
# Name: curvature.im
#
# Usage: curvature.im 
#
# Author: Yasin Dagasan
#
# Description:
# -----------
# This R function is used to compute different curvature types for a given image.
# The curvatures are computed as described in the following publication:
# -----------------------------------------------------------------------------
# Roberts, A. (2001). Curvature attributes and their application to 3D 
# interpreted horizons. First break, 19(2), 85-100.
# -----------------------------------------------------------------------------
#
# Input:  2D image in a matrix form
# Output: list of different types of curvature maps
###############################################################################

#Node numbering for the 3x3 window
#Z1    Z2    Z3
#Z4    Z5    Z6
#Z7    Z8    Z9

curvature.im<- function(image, gridsize, mask=F, mask.im=NULL){
  
  require(Matrix)  
  matdim<-dim(image)
  
  # Generation of curvature maps consisting of NA values
  mean.curv.map     <-matrix(NA,matdim[1], matdim[2])
  gaussian.curv.map <-matrix(NA,matdim[1], matdim[2])
  maximum.curv.map  <-matrix(NA,matdim[1], matdim[2])
  minimum.curv.map  <-matrix(NA,matdim[1], matdim[2])
  mostpos.curv.map  <-matrix(NA,matdim[1], matdim[2])
  mostneg.curv.map  <-matrix(NA,matdim[1], matdim[2])
  shape.index.map   <-matrix(NA,matdim[1], matdim[2])
  contour.curv.map  <-matrix(NA,matdim[1], matdim[2])
  dip.curv.map      <-matrix(NA,matdim[1], matdim[2])
  strike.curv.map   <-matrix(NA,matdim[1], matdim[2])
  curvedness.map    <-matrix(NA,matdim[1], matdim[2])
  
  # size of the kernel
  template.size<-c(3,3)
  
  # get the dimensions of the given image
  xdim = template.size[1]
  ydim = template.size[2]
  
  # construct the kernel
  lag = t(as.matrix(expand.grid( x = 0:(xdim - 1), y = 0:(ydim - 1))))
  n = dim(image)
  n <- c(n[1] - (xdim - 1), n[2] - (ydim - 1))
  pat = matrix(0, n[1] * n[2], dim(lag)[2])
  
  k = 1
  for (i2 in 1:n[1]) {
    for (i1 in 1:n[2])  {
      u = c(i2, i1)
      v = u + lag
      t = v+c(1,1)
        pat[k,] = image[(v[2,] - 1) * dim(image)[1] + v[1,]]
        
        if (!anyNA(pat[k,])) {
          z<-pat[k,]
          z1<-z[1]
          z2<-z[4]
          z3<-z[7]
          z4<-z[2]
          z5<-z[5]
          z6<-z[8]
          z7<-z[3]
          z8<-z[6]
          z9<-z[9]
          
          a <- (z1+z3+z4+z6+z7+z9)/12*(gridsize^2)-(z2+z5+z8)/6*(gridsize^2)
          b <- (z1+z2+z3+z7+z8+z9)/12*(gridsize^2)-(z4+z5+z6)/6*(gridsize^2)
          c <- (z3+z7-z1-z9)/4*(gridsize^2)
          d <- (z3+z5+z9-z1-z4-z7)/6*(gridsize^2)
          e <- (z1+z2+z3-z7-z7-z9)/6*(gridsize^2)
          f <- (2*(z2+z4+z6+z8)-(z1+z3+z7+z9)+5*z5)/9
          
          Km   <- (a*(1+e^2)-c*d*e+b*(1+d^2))/(1+d^2+e^2)^(3/2)
          Kg   <- (4*a*b-c^2)/(1+d^2+e^2)^2
          Kmax <- Km+sqrt(Km^2-Kg)
          Kmin <- Km-sqrt(Km^2-Kg)
          Kpos <- (a+b)+sqrt((a-b)^2+c^2)
          Kneg <- (a+b)-sqrt((a-b)^2+c^2)
          Si   <- (2/pi)*atan((Kmin+Kmax)/(Kmax-Kmin))
          Kd   <- 2*(a*e^2-c*d*e+b*d^2)/((d^2+e^2)*(1+d^2+e^2))^(1/2)
          Ks   <- 2*(a*e^2-c*d*e+b*d^2)/((d^2+e^2)*(1+d^2+e^2))^(1/2)
          Kc   <- 2*(a*e^2-c*d*e+b*d^2)/(1+d^2+e^2)^(3/2)
          Kn   <- sqrt((Kmax^2+Kmin^2)/2)
          
          mean.curv.map[t[1],t[2]]     <- Km
          gaussian.curv.map[t[1],t[2]] <- Kg
          maximum.curv.map[t[1],t[2]]  <- Kmax
          minimum.curv.map[t[1],t[2]]  <- Kmin
          mostpos.curv.map[t[1],t[2]]  <- Kpos
          mostneg.curv.map[t[1],t[2]]  <- Kmin
          shape.index.map[t[1],t[2]]   <- Si
          dip.curv.map[t[1],t[2]]      <- Kd
          strike.curv.map[t[1],t[2]]   <- Ks
          contour.curv.map[t[1],t[2]]  <- Kc
          curvedness.map[t[1],t[2]]    <- Kn
          
        }
        
        else {
          mean.curv.map[t[1],t[2]]     <- NA
          gaussian.curv.map[t[1],t[2]] <- NA
          maximum.curv.map[t[1],t[2]]  <- NA
          minimum.curv.map[t[1],t[2]]  <- NA
          mostpos.curv.map[t[1],t[2]]  <- NA
          mostneg.curv.map[t[1],t[2]]  <- NA
          shape.index.map[t[1],t[2]]   <- NA
          dip.curv.map[t[1],t[2]]      <- NA
          strike.curv.map[t[1],t[2]]   <- NA
          contour.curv.map[t[1],t[2]]  <- NA
          curvedness.map[t[1],t[2]]    <- NA
          
          
        }
        
        k = k + 1
        
    }
  }
  
  if (mask) {
    mean.curv.map[which(is.na(mask.im))]     <- NA
    gaussian.curv.map[which(is.na(mask.im))] <- NA
    maximum.curv.map[which(is.na(mask.im))]  <- NA
    minimum.curv.map[which(is.na(mask.im))]  <- NA
    mostpos.curv.map[which(is.na(mask.im))]  <- NA
    mostneg.curv.map[which(is.na(mask.im))]  <- NA
    shape.index.map[which(is.na(mask.im))]   <- NA
    dip.curv.map[which(is.na(mask.im))]      <- NA
    strike.curv.map[which(is.na(mask.im))]   <- NA
    contour.curv.map[which(is.na(mask.im))]  <- NA
    curvedness.map[which(is.na(mask.im))]    <- NA
  }
  
  
  result<-list(  mean.curve    = mean.curv.map,
                 gaussian.curv = gaussian.curv.map,
                 maximum.curv  = maximum.curv.map,
                 minimum.curv  = minimum.curv.map,
                 mostpos.curv  = mostpos.curv.map,
                 mostneg.curv  = mostneg.curv.map,
                 shape.index   = shape.index.map,
                 dip.curv      = dip.curv.map,
                 strike.curv   = strike.curv.map,
                 contour.curv  = contour.curv.map,
                 curvedness    = curvedness.map                 
                 )
  
  return(result)
  
}