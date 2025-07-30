Ispline<-function(x,order,knots){
  # M Spline function with order k=order+1. or I spline with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.
  
  k=order+1
  m=length(knots)
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots
  
  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }
  
  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }
  
  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }
  
  yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))
  
  if (order==1){
    for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          yy[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          yy[i-1,j]=0
        }
      }
    }
  }
  return(yy)
}

#files <- list.files(path = "F:/pro_left_trun",pattern = 'rtem.........txt|rtem........txt|rtem.......txt') 
#dataset <- list.files(path = "F:/pro_left_trun",pattern = 'data......txt|data........txt|data.......txt')

files <- list.files(path = "C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data",pattern = 'rtemr.....txt|rtemr....txt|rtemr...txt') 
dataset <- list.files(path = "C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data",pattern = 'datar.....txt|datar....txt|datar...txt')

setwd('C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data')
#number of observations per set
n <- 200
#the number of interior knots
ik <- 9
#Degree for Ispline
dgr <- 3
#number of parameters
K <- dgr + ik

#matrix to restore the result
result <- matrix(0,length(files),13)
#9 plots in one picture
par(mfrow=c(2,2))

for(j in files) {
  data_name <- dataset[which(files == j)]
  mydata <- read.table(paste('C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data',data_name,sep = "/"), header = TRUE, sep = "\t")
  rt <- read.table(j, header=TRUE)
  #real beata
  realbt <- c(round(mean(rt[,1])),round(mean(rt[,2])))
  #time points
  t <- seq(0,6,0.05)
  #the value of CDF for every t
  bslodds <- matrix(rep(NA,length(t)*500),500,length(t))
  
  for(i in 1:500){
    #defind data in each data set
    fst <- (i-1)*n+1
    lst <- i*n
    D <- mydata[fst:lst, ]
    
    #equal space knots
    #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
    mx <-  max(setdiff(c(D[,"R"],D[,"L"]),Inf)) + 0.01
    knots <- seq(0,mx,length.out = ik + 2)
    
    #calculate baseline odds
    gm <- rt[i,3:(2+K)]
    b <- Ispline(t,dgr,knots)
    yy <- as.matrix(gm)%*% b
    
    bslodds[i,] <- yy
  }
  
  odds <- apply(bslodds,2,mean)
  
  ############################need to change if baseline odds changed##############################
  bsl <- log(1+t) + t^3 + sin(t)
  #bsl <- log(1+t) + t^1.5
  
  #check baseline odds
  plot(t,odds,cex=0.5,ylab = "baseline odds")
  points(t,bsl,type = "l")
  
  #check baseline cdf
  plot(odds/(1+odds),bsl/(1+bsl),type = "p",col="blue",xlab = "estimated CDF", ylab = "true CDF")
  abline(a = 0,b  = 1,col = "black")
  plot(t,odds/(1+odds),ylab = "baseline CDF",ylim = c(0,1))
  points(t,bsl/(1+bsl),type = "l")
  
  mse <- apply((t(bslodds/(1+bslodds)) - bsl/(1+bsl))^2,1,mean)
  #check baseline odds
  plot(t,mse,cex=0.5,ylab = "mse",pch = 19,main = j)
  
  #check beta mse
  bt1_mse <- round(mean((rt[,1] - realbt[1])^2),4)
  bt2_mse <- round(mean((rt[,2] - realbt[2])^2),4)
  
  rt <- rt[complete.cases(rt),]
  
  #bias
  bias <- round(apply(rt,2,mean)[c(1,2)] - realbt,4)
  #check 95% ci
  ci_1 <- sum(rt[,17] < realbt[1] & rt[,18] > realbt[1])/nrow(rt)
  ci_2 <- sum(rt[,19] < realbt[2] & rt[,20] > realbt[2])/nrow(rt)
  #se
  esd <- apply(rt,2,mean)[c(15,16)]
  #mean running time
  meanrt <- round(mean(rt[,21]),4)
  
  result[which(files == j),] <- c(bias,esd,sd(rt[,1]),sd(rt[,2]),ci_1,ci_2,bt1_mse,bt2_mse,mean(mse),max(mse),meanrt)
  
  A <- cbind(files,result)
}

library(xtable)
newobject2 <- as.data.frame(A)
xtable(newobject2, type = "latex", file = "filename2.tex")

#write.csv(newobject2,file = 'F:/pro_left_trun/result.csv')
write.csv(newobject2,file = 'C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data/resultr.csv')


