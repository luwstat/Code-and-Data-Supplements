
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
dataset <- list.files(path = "C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data",pattern = 'datar.........txt|datar........txt|datar.......txt')
setwd('C:/Users/wangl/Dropbox/projectLeftTruncation/simulation/data')
#set up a vector of result names
resultsets <- gsub("data","rtema",dataset)

#number of observations per set
n <- 200
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
    I[,D[,"d_2"] == 1] <- Ir[,D[,"d_2"] == 1] - Il[,D[,"d_2"] == 1]
    
    #set an initial value
    bt <- rep(1,P)
    gama <- rep(1,K)
    df <- rep(1,P)
    ite <- 0
    Ez <- rep(0,n)
    Ew <- rep(0,n)
    Epsi <- rep(0,n)
    Ephi <- rep(0,n)
    EU <- matrix(0,n,K)
    
    while( t(df)%*%df > 1e-8 & ite < 20000){
      #exp(x_i*beta) is n*1 vector
      exb <- exp(iv%*%bt)
      #Lambda_0(R) and Lambda_0(L), n*1 vector
      Bsr <- t(Ir)%*%gama
      Bsl <- t(Il)%*%gama
      #Lambda_0(.)exp(x*beta) + 1, n*1 vector
      Bre <- Bsr*exb + 1
      Ble <- Bsl*exb + 1
      
      #Ez does not equal to zero only when d_1=1
      Ez[d_1 == 1] <- Bre[d_1 ==1]
      #Ew does not equal to zero only when d_2=1
      Ew[d_2 == 1] <- Bre[d_2 ==1]/Ble[d_2 ==1]
      #Ephi and Epsi are n*1 vectors
      Ephi[d_0 == 1 | d_3 ==1] <- 1/Ble[d_0 == 1 | d_3 ==1]
      Ephi[d_1 == 1] <- (Bre[d_1 == 1] + 1)/Bre[d_1 == 1]
      Ephi[d_2 == 1] <- (Bre + Ble)[d_2 == 1]/(Ble*Bre)[d_2 == 1]
      #Epsi does not equal to zero only when delta_0 = 1
      Epsi[d_0 == 1] <- 1/Ble[d_0 ==1]   
      #EZ and EW are n*K matrix
      EZ <- t(Ir*gama)*as.vector(Ez/Bsr)
      EW <- t(I*gama)*as.vector(Ew/(Bsr - Bsl)^d_2)
      #EU is a n*K matrix, the ^del is used to get rid of the 0 denomator in right censored data
      EU[d_0 == 1,] <- t(Ml[,d_0 == 1]*gama)/as.vector(t(Ml[,d_0 == 1])%*%gama)   
      
      #set the equations needed to be solved
      #first define the constant terms
      #A is a K*1 vector
      A <- apply( as.matrix(EU[d_0 == 1,]),2,sum) + apply( EZ[d_1 == 1,],2,sum) + apply( EW[d_2 == 1,],2,sum)
      #B is a n*K matrix
      B <- t(Il)*Epsi + ( t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2) )*Ephi
      one <- rep(1,K)
      btf <- function(x){
        y <- numeric(length(bt))
        for(h in 1:length(bt)){
          y[h] <- (d_0 + Ez + Ew)%*%iv[,h] - (one%*%((t(B)*A)/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
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
    #Lambda_0(.)exp(x*beta) + 1, n*1 vector
    Bre <- Bsr*exb + 1
    Ble <- Bsl*exb + 1
    
    #Ez not zero when d_1=1
    Ez[d_1 == 1] <- Bre[d_1 ==1]
    #Ew not zero when d_2=1
    Ew[d_2 == 1] <- Bre[d_2 ==1]/Ble[d_2 ==1]
    #Ephi and epsi are n*1 vectors
    Ephi[d_0 == 1 | d_3 ==1] <- 1/Ble[d_0 == 1 | d_3 ==1]
    Ephi[d_1 == 1] <- (Bre[d_1 == 1] + 1)/Bre[d_1 == 1]
    Ephi[d_2 == 1] <- (Bre + Ble)[d_2 == 1]/(Ble*Bre)[d_2 == 1]
    #Epsi only when delta_0 = 1
    Epsi[d_0 == 1] <- 1/Ble[d_0 ==1]   
    #EZ and EW are n*K matrix
    EZ <- t(Ir*gama)*as.vector(Ez/Bsr)
    EW <- t(I*gama)*as.vector(Ew/(Bsr - Bsl)^d_2)
    #EU is a n*K matrix, the ^del is used to get rid of the 0 denomator
    EU[d_0 == 1,] <- t(Ml[,d_0 == 1]*gama)/as.vector(t(Ml[,d_0 == 1])%*%gama)
    
    #########variances#############
    #variance phi is n*1 vector
    Vphi <- rep(0,n)
    Vphi[d_0 == 1 | d_3 ==1] <- Ephi[d_0 == 1 | d_3 ==1]^2
    Vphi[d_1 == 1] <- Bre[d_1 == 1]^(-2) + 1
    Vphi[d_2 == 1] <- Ble[d_2 == 1]^(-2) + Bre[d_2 == 1]^(-2)
    #variance psi is n*1 vector
    Vpsi <- rep(0,n)
    Vpsi[d_0 == 1]  <- Epsi[d_0 == 1]^2
    #Varicnce of z not zero when d_1=1 (vector of Var(z_i))
    Vz <- rep(0,n)
    Vz[d_1 == 1] <- Bre[d_1 == 1]*(Bre[d_1 == 1] - 1)
    #Varicnce of w not zero when d_2=1 (vector of Var(w_i))
    Vw <- rep(0,n)
    Vw[d_2 == 1] <- Ew[d_2 == 1]*(Ew[d_2 == 1] - 1  )
    #Variance of Z and W are n*K matrix
    VZ <- EZ*(1 - t(Ir*gama)/as.vector(Bsr)) + Vz*(t(Ir*gama)/as.vector(Bsr))^2
    VW <- EW*(1 - t(I*gama)/as.vector(Bsr - Bsl)^d_2) + Vw*(t(I*gama)/as.vector(Bsr - Bsl)^d_2)^2
    #VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
    VU <- EU*(1- EU)
    
    ####part 1 Q########################################
    Q <- matrix(0,(K+P),(K+P))
    #same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
    A <- apply( as.matrix(EU[d_0 == 1,]),2,sum) + apply( EZ[d_1 == 1,],2,sum) + apply( EW[d_2 == 1,],2,sum)
    B <- t(Il)*Epsi + ( t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2) )*Ephi
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
    D_i <- Bsl*(d_0 + d_3) + Bsr*(d_1 + d_2)
    D_il <- t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2)
    
    #cov(l_c/beta,l_c/beta )
    e_covc <- Vpsi*(Bsl*exb*d_0)^2 + Vphi*(D_i*exb)^2 + Vz*d_1 + Vw*d_2 - 2*(Bre - 1)^2*d_1 - 2*((Bre - Ble)/Ble^2)*(Bre - 1)*d_2
    covc <-  rep(e_covc, each = P^2)
    vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P) 
    
    #coefficients for cov(l_c/beta,l_c/gama )
    co1 <- t(Il)*as.vector(Bsl*exb^2)*Vpsi*d_0 + D_il*as.vector(D_i*exb^2)*Vphi
    co2 <- t(Ir)*as.vector(d_1*exb*Bre) + t(I)*as.vector( Ew*d_2*exb/Ble )
    co3 <- t(Ir)*as.vector(d_1*Bsr*exb^2) + t(Ir)*as.vector(d_2*(Bsr - Bsl)*exb^2*Ble^(-2)) + t(Ir)*as.vector(Bsr*exb^2*d_1) + t(I)*as.vector(d_2*Bsr*exb^2/Ble^2)
    co_t <- co1 + co2 - co3
    
    vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
    vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])
    #coefficients for cov(l_c/gama, l_c/gama) on diagnal 
    dg1 <- t(VU*d_0 + VW*d_2 + VZ*d_1)/gama^2
    dg2 <- ( t(Il^2)*Vpsi*d_0 + Vphi*D_il^2 )*as.vector(exb^2)
    dg3 <- t(Ir)*t(I)*d_2*as.vector(exb^2/(Ble^2)) + t(Ir^2)*d_1*as.vector(exb^2)
    
    diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2 -2*dg3,2,sum)
    #coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal 
    for(l in 1:K){
      for(m in 1:K){
        if(l != m){
          part_1 <- -(Ml[l,]*Ml[m,]*d_0/as.vector(t(Ml)%*%gama)^(2*d_0)) + d_1*Ir[l,]*Ir[m,]*(Vz-Ez)/Bsr^2 + d_2*I[l,]*I[m,]*(Vw - Ew)/(Bsr - Bsl)^(2*d_2)
          part_2 <- (Il[l,]*Il[m,]*d_0*Vpsi + D_il[,l]*D_il[,m]*Vphi)*exb^2
          part_3 <- d_1*Ir[m,]*Ir[l,]*exb^2 + d_2*Ir[m,]*I[l,]*exb^2/Ble^2
          part_4 <- d_1*Ir[l,]*Ir[m,]*exb^2 + d_2*Ir[l,]*I[m,]*exb^2/Ble^2
          vc[P+l,P+m] <- sum( part_1 + part_2 - part_3 - part_4 )
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
    
    
    lb <- bt - 1.96*se_theta
    ub <- bt + 1.96*se_theta
    ci <- as.vector(rbind(lb,ub))
    
    fittime <- (proc.time() - ptm)["elapsed"]
    
    result[i,] <- c(bt,gamanew,se_theta,ci,fittime)
    
    print(round(c(i,ite,bt)))
    
    #stop the time
    print(proc.time() - ptm)
    
  }
  result_file <- resultsets[which(dataset == j)]
  write.table(result, result_file, sep="\t",row.names = FALSE)
  
}

