
files <- list.files(path = "F:/pro_left_trun",pattern = 'rticenr_200....txt|rticenr_200.....txt|rticenr_200...txt')
dataset <- list.files(path = "F:/pro_left_trun",pattern = 'datar_200.....txt|datar_200....txt|datar_200...txt')
setwd("F:/pro_left_trun")
#number of observations per set
n <- 200

#matrix to restore the result
result <- matrix(0,length(files),13)
#9 plots in one picture
par(mfrow=c(2,2))

t <- seq(0,6,0.05)
odd <- log(1+t) + t^1.5
#odd <- log(1+t) + t^3 +sin(t)
yy <- 1/(1+odd)

for(j in files){
  data_name <- dataset[which(files == j)]
  mydata <- read.table(paste('F:/pro_left_trun',data_name,sep = "/"), header = TRUE, sep = "\t")
  rt1 <- read.table(j, header=TRUE)
  rt <- rt1[complete.cases(rt1),]
  
  #real beata
  realbt <- c(round(mean(rt[,1])),round(mean(rt[,2])))
  #check beta mse
  bt1_mse <- round(mean((rt[,1] - realbt[1])^2),4)
  bt2_mse <- round(mean((rt[,2] - realbt[2])^2),4)
  
  #bias
  bias <- round(apply(rt,2,mean)[c(1,2)] - realbt,4)
  #check 95% ci
  ci_1 <- sum(rt[,5] < realbt[1] & rt[,6] > realbt[1])/nrow(rt)
  ci_2 <- sum(rt[,7] < realbt[2] & rt[,8] > realbt[2])/nrow(rt)
  #se
  esd <- apply(rt,2,mean)[c(3,4)]
  
  #the value of S(t) for every t
  S <- rt[,seq(9,129)]
  
  mse <- apply((t(S) - yy)^2,1,mean)
  
  #check baseline odds
  plot(t,mse,cex=0.5,ylab = "mse",pch = 19,main = j)
  
  #mean running time
  meanrt <- round(mean(rt[,130]),4)
  
  result[which(files == j),] <- c(bias,esd,sd(rt[,1]),sd(rt[,2]),ci_1,ci_2,bt1_mse,bt2_mse,mean(mse),max(mse),meanrt)
  
  A <- cbind(files,result)
  
}

library(xtable)
newobject2 <- as.data.frame(A)
xtable(newobject2, type = "latex", file = "filename2.tex")

write.csv(newobject2,file = 'F:/pro_left_trun/rticenr.csv')

