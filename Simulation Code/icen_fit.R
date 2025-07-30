library(icenReg)

#read data automatically
dataset <- list.files(path = "F:/pro_left_trun",pattern = 'datar_200.....txt|datar_200....txt|datar_200...txt')

setwd('F:/pro_left_trun')
#set up a vector of result names
resultsets <- gsub("data","rticen",dataset)

t <- seq(0,6,0.05)

for (j in dataset) {
  #read data in each datasets
  mydata <- read.table(j, header = TRUE, sep = "\t")
  result <- NULL
  #number of observations per set
  n <- nrow(mydata)/500
  
  for(i in 1:500){
    #start the clock
    ptm <- proc.time()
    
    #choose data in each data set
    fst <- (i-1)*n+1 
    lst <- i*n
    D <- mydata[fst:lst, ]
    
    left <- D[,"L"]
    right <- D[,"R"]
    fit_po <- ic_sp(cbind(left,right)~x_1+x_2,model = 'po',bs_samples = 100,data = D)
    
    bt <- -fit_po$coefficients
    se_theta <- sqrt(diag(fit_po$var))
    
    lb <- bt - 1.96*se_theta
    ub <- bt + 1.96*se_theta
    ci <- as.vector(rbind(lb,ub))
    
    fittime = (proc.time() - ptm)["elapsed"]
    bs <- data.frame(x_1 = 0, x_2 = 0)
    A <- getSCurves(fit_po,newdata = bs)
    st <- c(1,A$S_curves$`1`[findInterval(t,A$Tbull_ints[,1])])
    
    #record the solutions
    result <- rbind(result,c(bt,se_theta,ci,st,fittime))
    print(c(i,bt))
    
  }
  
  result_file <- resultsets[which(dataset == j)]
  write.table(result,paste("F:/pro_left_trun",result_file,sep = "/") , sep="\t",row.names = FALSE)
}
