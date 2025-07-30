#install.packages("devtools") #install this pacakge if you haven’t installed it yet.
#library(devtools)
# Install the proposed R pacakge regPOspline if you haven't installed it yet.
#remotes::install_github("luwstat/regPOspline", dependencies = TRUE)
library(regPOspline)

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
#read data
load("H:/My Drive/Dropbox/projectLeftTruncation/realdata/finalDataInfection3LT.rdata")

#number of observations per set
n <- nrow(dataInfection3LTbirth)
# the number of lifetime partners reported at the time of enrollment 
x_1 = dataInfection3LTbirth[,"nptrlife0"]
# 1: African American 0: otherwise
x_2 = as.numeric(dataInfection3LTbirth[,"race"] == 1)
C = dataInfection3LTbirth[,"agefrstsx"]

L <- dataInfection3LTbirth[,c("L_TV")]
R <- replace(dataInfection3LTbirth[,c("R_TV")],is.na(dataInfection3LTbirth[,c("R_TV")]), Inf)

d_1 <- rep(0, n)
d_2 <- rep(0, n)
d_3 <- rep(0, n)
d_0 <- as.numeric(L == R)
for (i in which(d_0 == 0)) {
  if (L[i] == C[i]) {
    d_1[i] <- 1
  }
  else if (R[i] == Inf) {
    d_3[i] <- 1
  }
  else {
    d_2[i] <- 1
  }
}

Data <- cbind(L,R,C,d_0,d_1,d_2,d_3,x_1,x_2)

set.seed(10)  # for reproducibility

# Create a random index
train_index <- sample(1:nrow(Data), 0.7 * nrow(Data))  # 70% for training

# Split the data
D <- Data[train_index, ]   # D: training data
Dt  <- Data[-train_index, ]   # Dt: testing data
#ivt is the X's in testing data
ivt <- as.matrix(Dt[,c("x_1","x_2")])
#the number of interior knots
ik <- 1
#Degree for Ispline
dgr <- 2
#number of parameters
K <- dgr + ik
mx <-  max(setdiff(c(Data[,"R"],Data[,"L"]),Inf)) + 0.01
knots <- seq(0,mx,length.out = ik + 2)


# Fit the model using training data D
result <- po_fit(
  L = D[, "L"],
  R = D[, "R"],
  truncation = F,
  x = D[, c("x_1","x_2")],
  initial_reg = runif(2,-2,2),
  order = dgr,
  equal_space = F,
  myknots = knots,
  initial_spline = rep(1, K)
)

# Save estimated beta in bt and estimated gamma in gama
bt = result$coefficient_est[,"bt"]
gama = result$spline_coef

# Calculate IBS based on the testing data
# Define time grids to approximate IBS
delta_t = 0.05
time_grid <- seq(0,mx,delta_t)
b <- Ispline(time_grid,dgr,knots)
# Get S(t|x) for each individual, so each individual has length(time_grid) survival probs
# Saving this in a matrix
odds <- exp(ivt%*%bt)%*%(t(as.matrix(gama)) %*% b)
sur <- 1/(1 + odds)

po_surv <- function(bdr,xi){
  Ibdr <- Ispline(bdr,dgr,knots)
  bdrodds <- exp(t(xi)%*%bt)%*%(t(as.matrix(gama)) %*% Ibdr)
  return(1/(1 + bdrodds))
}
# Function to compute IBS for a single subject i
compute_IBS_i <- function(i) {
  Li <- Dt[i,"L"]
  Ri <- Dt[i,"R"]
  Xi <- as.matrix(Dt[i,c("x_1","x_2")])
  sq_errors <- numeric(length(time_grid))
  
  for (j in seq_along(time_grid)) {
    t <- time_grid[j]
    S_hat_t <- sur[i,j]
    
    # Determine I_hat(T_i > t)
    if (t < Li) {
      I_hat <- 1
    } else if (t >= Ri) {
      I_hat <- 0
    } else {
      # Avoid divide-by-zero if S(Li) == S(Ri)
      denom <- po_surv(Li,Xi) - po_surv(Ri,Xi)
      if (abs(denom) < 1e-8) {
        I_hat <- 0.5  # or some neutral value
      } else {
        I_hat <- (S_hat_t - po_surv(Ri,Xi)) / denom
        I_hat <- min(max(I_hat, 0), 1)  # clip to [0, 1]
      }
    }
    
    # Squared error
    sq_errors[j] <- (I_hat - S_hat_t)^2
  }
  
  # Approximate integral for subject i using average (Riemann)
  mean(sq_errors)
}

# Compute IBS across all subjects
IBS_all <- sapply(1:nrow(Dt), compute_IBS_i)
IBS <- mean(IBS_all)
cat("Estimated IBS:", IBS, "\n")








