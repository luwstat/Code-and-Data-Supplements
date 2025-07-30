#install.packages("devtools") #install this pacakge if you haven’t installed it yet.
library(devtools)
# Install the proposed R pacakge regPOspline if you haven’t installed it yet.
#remotes::install_github("luwstat/regPO", dependencies = TRUE)
library(regPO)
####################Read data set and change setting here####################
#data cleaning
my.data <- read.table(file = "F:/project1/revise_realdata/C.csv", sep = "," ,header = TRUE)

#mutant type
my.data <- my.data[my.data[,"KRAS_C"] == 1,]
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

sp_fit(L=D[,"L"], R=D[,"R"], x=D[,"TRT_C"], order =2, nknot = 11 )


#data cleaning
my.data <- read.table(file = "F:/project1/revise_realdata/C.csv", sep = "," ,header = TRUE)
#wild type
my.data <- my.data[my.data[,"KRAS_C"] == 0,]

L <- rep(NA,n)
R <- rep(NA,n)
L[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0,1]
R[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0,1]

L[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1,"L"]
R[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1,"R"]
R[my.data[,"y"] == 2] <- Inf

TRT_C <- my.data[,c("TRT_C")]
D <- cbind(L,R,TRT_C)

sp_fit(L=D[,"L"], R=D[,"R"], x=D[,"TRT_C"], order =2, nknot = 3 )