
###################Functions will be used in the EM algorithm###############
library(nleqslv)
library(Matrix)
#M and I-Spline function
Mspline<-function(x,order,knots){
  # get M spline and I spline matrix with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.
  
  
  ### get Mspline bases ###
  k1=order
  m=length(knots)
  n1=m-2+k1 # number of parameters
  t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # newknots
  
  tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
  for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
  }
  
  if (order==1){
    mbases=tem1
  }else{
    mbases=tem1
    for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
        tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
    }
  }
  
  return(mbases)
}

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


####################Read data set and change setting here####################
#read data automatically
#dataset <- list.files(path = "F:/pro_left_trun",pattern = 'dataold....txt|dataold.....txt|dataold...txt')

dataset <- list.files(path = "C:/Users/LW HomeOffice/Dropbox/projectLeftTruncation/simulation/data",pattern = 'data.........txt|data........txt|data.......txt')

setwd('C:/Users/LW HomeOffice/Dropbox/projectLeftTruncation/simulation/data')
#set up a vector of result names
# rtem restult when initial values are fixed to 1
#resultsets <- gsub("data","rtem",dataset)
# rtemrd restult when initial values are random
resultsets <- gsub("data","rtem_rdm",dataset)

#the number of interior knots
ik <- 9
#Degree for Ispline
dgr <- 3
#number of parameters
K <- dgr + ik
##################EM algorithm##############################
for (j in dataset) {
  #read data in each datasets 
  mydata <- read.table(j, header = TRUE, sep = "\t")
  #restore result for each datasets in a matrix
  result <- matrix(rep(NA,(9+K)*500),500,(9+K))
  #number of observations per set
  n <- nrow(mydata)/500
  for(i in 1:500){
    #start the clock
    ptm <- proc.time()
    
    #choose data in each data set
    fst <- (i-1)*n+1
    lst <- i*n
    D <- mydata[fst:lst, ]
    #iv is the X's
    iv <- as.matrix(D[,c("x_1","x_2")])
    P <- ncol(iv)
    #del is delta
    d_0 <- D[,"d_0"]
    d_1 <- D[,"d_1"]
    d_2 <- D[,"d_2"]
    d_3 <- D[,"d_3"]
    
    #equal space knots
    #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
    mx <-  max(setdiff(c(D[,"R"],D[,"L"]),Inf)) + 0.01
    knots <- seq(0,mx,length.out = ik + 2)
    
    Ml <- Mspline(D[,"L"],dgr,knots)
    Mr <- Mspline(D[,"R"],dgr,knots)
    Il <- Ispline(D[,"L"],dgr,knots)
    Ir <- Ispline(D[,"R"],dgr,knots)
    #rate for interval censored data
    I <- matrix(0,K,n)
    I[,d_1 == 1 | d_2 == 1] <- Ir[,d_1 + d_2 == 1] - Il[,d_1 == 1 | d_2 == 1]
    
    #left truncation time spline
    Ic <- Ispline(D[,"C"],dgr,knots)
    
    #set an initial value
    #bt <- rep(0,P)
    bt <- runif(P,-2,2)
    gama <- rep(1,K)
    
    ################################
    df <- rep(1,P)
    ite <- 0
    
    while( max(abs(df)) > 1e-5 & ite < 20000){
      #exp(x_i*beta) is n*1 vector
      exb <- exp(iv%*%bt)
      #Lambda_0(R) and Lambda_0(L), n*1 vector
      Bsr <- t(Ir)%*%gama
      Bsl <- t(Il)%*%gama
      
      #Ephi and Epsi are n*1 vectors
      Ephi <- as.vector(1/(Bsl*exb + 1))
      Epsi <- as.vector(1/(Bsr*exb + 1))
      Epsi[d_3 == 1] = rep(0,sum(d_3))
 
      #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
      EU <- t(I*gama)/as.vector(t(I)%*%gama)   
      EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)   
      EU[d_3 == 1,] <- rep(0,K)
      
      #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
      EV_1 <- 1/(1 + as.vector(t(Ic)%*%gama)*exb)
      EV_2 <- (t(Ic*gama)*as.vector(exb))/as.vector(1 + as.vector(t(Ic)%*%gama)*exb)
      EV <- cbind(EV_1,EV_2)
      
      #The equations needed to be solved, x in the function to be solved is the beta
      #A is a K*1 vector
      A <- colSums(EV_2) + colSums(EU)
      #B is a n*K matrix
      B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi
      
      #one <- rep(1,K)
      #btf <- function(x){
      #  y <- numeric(length(bt))
      #  for(h in 1:length(bt)){
      #    y[h] <- (d_0 + d_1 + d_2 + rowSums(EV_2))%*%iv[,h] - (one%*%((t(B)*A)/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
      # }
      #  y
      #}
      btf <- function(x){
        y <- numeric(length(bt))
        for(h in 1:length(bt)){
         y[h] <- (d_0 + d_1 + d_2 + rowSums(EV_2))%*%iv[,h] - t(B%*%(A/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
        }
        y
      }
      btstart <- rep(0,length(bt))
      #solve and get the updated bt
      sol <- nleqslv(btstart,btf,method="Newton")
      btnew <- sol$x
      
      gamanew <-  A/( as.vector(t(B)%*%exp(iv%*%btnew)) )
      
      df <- btnew - bt
      bt <- as.vector(btnew)
      gama <- as.vector(gamanew)
      ite <- ite + 1
      
    }
    
    
    ################calculate covariance matrix##################
    #First calculate all the expectationes and variances given the theta_hat
    ######Expectations######
    #x_i*beta is n*1 vector
    xb <- iv%*%bt
    #exp(x_i*beta) is n*1 vector
    exb <- exp(xb)
    #Lambda_0(R) and Lambda_0(L), n*1 vector
    Bsr <- t(Ir)%*%gama
    Bsl <- t(Il)%*%gama

    #Ephi and Epsi are n*1 vectors
    Ephi <- as.vector(1/(Bsl*exb + 1))
    Epsi <- as.vector(1/(Bsr*exb + 1))
    Epsi[d_3 == 1] = rep(0,sum(d_3))
    
    #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
    EU <- t(I*gama)/as.vector(t(I)%*%gama)   
    EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)   
    EU[d_3 == 1,] <- rep(0,K)
    
    #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
    EV_1 <- as.vector( 1/(1 + as.vector(t(Ic)%*%gama)*exb) )
    EV_2 <- (t(Ic*gama)*as.vector(exb))/as.vector(1 + as.vector(t(Ic)%*%gama)*exb)

    
    #########variances#############
    #variance phi is n*1 vector
    Vphi <- as.vector(1/(Bsl*exb + 1))^2
    Vpsi <- as.vector(1/(Bsr*exb + 1))^2
    Vpsi[d_3 == 1] = rep(0,sum(d_3))
    
    #VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
    VU <- EU*(1- EU)
    VV <- EV_2*(1- EV_2)
    
    ####part 1 Q########################################
    Q <- matrix(0,(K+P),(K+P))
    #same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
    A <- colSums(EV_2) + colSums(EU)
    B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi
    
    #elements in xx^t as a vector, four entries as a group
    e <- as.vector(t(cbind(iv,iv))) * rep(as.vector(t(iv)),each = P)
    #coefficients for each matrix, repeat each elememt p times so it match the vector above
    co <-  rep(B%*%gama*exb, each = P^2)
    #fill in the corresponding place in Q
    Q[seq(1,P),seq(1,P)] <- - matrix(apply(matrix(e*co,P^2,n),1,sum),P,P)
    
    Q[seq(P+1,K+P),seq(1,P)] <- - t(B)%*%(iv*as.vector(exb))
    
    Q[seq(1,P),seq(P+1,P+K)] <- t( Q[seq(P+1,K+P),seq(1,P)])
    
    diag(Q[seq(1+P,P+K),seq(1+P,P+K)]) <- - gama^(-2)*A
    
    ####part 2 VC ###########################################
    vc <- matrix(0,(K+P),(K+P))
    
    #cov(l_c/beta,l_c/beta )
    e_covc <- Vpsi*(Bsr*exb)^2 + Vphi*(Bsl*exb)^2 + EV_1*(1 - EV_1)
    covc <-  rep(e_covc, each = P^2)
    vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P) 
    
    #coefficients for cov(l_c/beta,l_c/gama )
    co_t <- t(Ir)*as.vector(Bsr*exb^2)*Vpsi + t(Il)*as.vector(Bsl*exb^2)*Vphi + t(t(EV_1*EV_2)/gama)
    vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
    vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])
    
    #coefficients for cov(l_c/gama, l_c/gama) on diagnal 
    dg1 <- t(VU + VV)/gama^2
    dg2 <- ( t(Ir^2)*Vpsi + Vphi*t(Il^2) )*as.vector(exb^2)
    
    diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2,2,sum)
    #coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal 
    for(l in 1:K){
      for(m in 1:K){
        if(l != m){
          part_1 <- -as.vector(EV_2[,l]*EV_2[,m] + EU[,l]*EU[,m])/as.vector(gama[l]*gama[m])
          part_2 <- (Ir[l,]*Ir[m,]*Vpsi + Il[l,]*Il[m,]*Vphi)*exb^2
          vc[P+l,P+m] <- sum( part_1 + part_2 )
        }
      } 
    }
    
    v <- -(Q + vc)
    
    ####part3
    tol <- 1e-7
    if(sum(gama < tol) == 0){
      vv <- v
    } else {
      rg=1:K
      index=rg[gama < tol]+P 
      vv=v[-index,-index]
    }
    
    if( rcond(vv) > .Machine$double.eps ){
      se_theta = sqrt(diag(solve(vv))[1:P])
    } else {
      se_theta = sqrt(diag(solve(vv + diag(1e-4,nrow(vv),ncol(vv))))[1:P])
    }
    
    
    ###############CI###########################
    lb <- bt - 1.96*se_theta
    ub <- bt + 1.96*se_theta
    ci <- as.vector(rbind(lb,ub))
    
    fittime <- (proc.time() - ptm)["elapsed"]
    
    result[i,] <- c(bt,gamanew,se_theta,ci,fittime)
    
    print(c(bt,i,ite))
    
    #stop the time
    print(proc.time() - ptm)
    
  }
  result_file <- resultsets[which(dataset == j)]
  write.table(result, result_file, sep="\t",row.names = FALSE)
  
}

