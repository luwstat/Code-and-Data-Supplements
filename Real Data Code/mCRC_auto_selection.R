####################Read data set and change setting here####################
#data cleaning
my.data <- read.table(file = "F:/project1/revise_realdata/C.csv", sep = "," ,header = TRUE)
my.data <- my.data[my.data[,"KRAS_C"] == 0,]
#my.data <- my.data[my.data[,"KRAS_C"] == 1,]
n <- nrow(my.data)

L <- rep(NA,n)
R <- rep(NA,n)
L[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0,1]
R[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0,1]

L[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1,"L"]
R[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1,"R"]
R[my.data[,"y"] == 2] <- Inf

TRT_C <- my.data[,c("TRT_C")]
D <- cbind(L,R,TRT_C)

result <- NULL
#############################################
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

for (dgr in 2:3) {
  for (ik in 1:15) {
    #number of parameters (length of gamma)
K <- dgr + ik
ptm <- proc.time()
mx <-  max(setdiff(c(D[,"R"],D[,"L"]),Inf)) + 0.01
knots <- seq(0,mx,length.out = ik + 2)
#rt <- sp_fit(L = D$L,R = D$R,x = D[,c(3,4)],order = 3,nknot = 3 )

iv <- as.matrix(D[,3])
P <- ncol(iv)
n <- nrow(iv)
d_1 <- rep(0, n)
d_2 <- rep(0, n)
d_3 <- rep(0, n)
d_0 <- as.numeric(L == R)
for (i in which(d_0 == 0)) {
  if (L[i] == 0) {
    d_1[i] <- 1
  }
  else if (R[i] == Inf) {
    d_3[i] <- 1
  }
  else {
    d_2[i] <- 1
  }
}
Ml <- Mspline(L, dgr, knots)
Mr <- Mspline(R, dgr, knots)
Il <- Ispline(L, dgr, knots)
Ir <- Ispline(R, dgr, knots)
#rate for interval censored data
I <- matrix(0,K,n)
I[,d_1 == 1 | d_2 == 1] <- Ir[,d_1 + d_2 == 1] - Il[,d_1 == 1 | d_2 == 1]
bt <- rep(1, P)
gama <- rep(1, K)
df <- rep(1, P)
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
  
  #The equations needed to be solved, x in the function to be solved is the beta
  #A is a K*1 vector
  A <- colSums(EU)
  #B is a n*K matrix
  B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi
  
  btf <- function(x){
    y <- numeric(length(bt))
    for(h in 1:length(bt)){
      y[h] <- (d_0 + d_1 + d_2)%*%iv[,h] - t(B%*%(A/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
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
#First calculate all the expectations and variances given the theta_hat
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

#########variances#############
#variance phi is n*1 vector
Vphi <- as.vector(1/(Bsl*exb + 1))^2
Vpsi <- as.vector(1/(Bsr*exb + 1))^2
Vpsi[d_3 == 1] = rep(0,sum(d_3))

#VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
VU <- EU*(1- EU)

####part 1 Q########################################
Q <- matrix(0,(K+P),(K+P))
#same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
A <- colSums(EU)
B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi

#elements in xx^t as a vector, four entries as a group
#e <- as.vector(t(iv)) * rep(as.vector(t(iv)),each = P)
e <- as.vector(t(do.call(cbind, replicate(P, iv, simplify = F)))) * rep(as.vector(t(iv)), each = P)
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
e_covc <- Vpsi*(Bsr*exb)^2 + Vphi*(Bsl*exb)^2
covc <-  rep(e_covc, each = P^2)
vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P) 

#coefficients for cov(l_c/beta,l_c/gama )
co_t <- t(Ir)*as.vector(Bsr*exb^2)*Vpsi + t(Il)*as.vector(Bsl*exb^2)*Vphi
vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])

#coefficients for cov(l_c/gama, l_c/gama) on diagnal 
dg1 <- t(VU)/gama^2
dg2 <- ( t(Ir^2)*Vpsi + Vphi*t(Il^2) )*as.vector(exb^2)

diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2,2,sum)
#coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal 
for(l in 1:K){
  for(m in 1:K){
    if(l != m){
      part_1 <- -as.vector(EU[,l]*EU[,m])/as.vector(gama[l]*gama[m])
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


llhd <- sum(log(((t(Ml) %*% gama) * exp(iv %*% bt)/(t(Il) %*% 
                                                      gama * exp(iv %*% bt) + 1)^2)^d_0 * (1 - 1/(t(Ir) %*% 
                                                                                                    gama * exp(iv %*% bt) + 1))^(d_1) * (1/(t(Il) %*% gama * 
                                                                                                                                              exp(iv %*% bt) + 1) - 1/(t(Ir) %*% gama * exp(iv %*% 
                                                                                                                                                                                              bt) + 1))^(d_2) * (1/(t(Il) %*% gama * exp(iv %*% bt) + 
                                                                                                                                                                                                                      1))^(d_3)))
AIC <- 2 * (P + K) - 2 * llhd
BIC <- (P+K)*log(n) - 2*llhd

fittime <- (proc.time() - ptm)["elapsed"]
result <- rbind(result,c(dgr,ik,bt,se_theta,as.vector(ci),AIC,BIC,fittime))
  }
  
}

write.csv(result, "F:/result_KRAS0.csv",row.names = FALSE)
