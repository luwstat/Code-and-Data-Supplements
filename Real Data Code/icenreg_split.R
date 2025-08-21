library(icenReg)

## helper to pretty-print a system.time object
print_timing <- function(tag, tm) {
  cat(sprintf("\n[%s] timing — elapsed: %.2f s | user: %.2f s | sys: %.2f s\n",
              tag, tm["elapsed"], tm["user.self"], tm["sys.self"]))
}

t_total_start <- proc.time()

#################### Read data set ####################
dat0 <- read.table(file = "F:/project1/revise_realdata/C.csv", sep = "," , header = TRUE)

#################### Mutant type (KRAS_C == 1) ####################
my.data <- dat0[dat0[,"KRAS_C"] == 1, ]
n <- nrow(my.data)

L <- R <- rep(NA_real_, n)
L[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0, 1]
R[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0, 1]
L[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1, "L"]
R[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1, "R"]
R[my.data[,"y"] == 2] <- Inf

TRT_C <- my.data[, "TRT_C"]
D <- cbind(L, R, TRT_C)

tm_mut <- system.time({
  fit0 <- ic_sp(cbind(D[,"L"], D[,"R"]) ~ TRT_C,
                model = 'po',
                data = as.data.frame(D),
                bs_samples = 100)
  s0 <- summary(fit0)
})
print(s0)
print_timing("KRAS mutant (KRAS_C==1)", tm_mut)

#################### Wild type (KRAS_C == 0) ####################
my.data <- dat0[dat0[,"KRAS_C"] == 0, ]
n <- nrow(my.data)

L <- R <- rep(NA_real_, n)
L[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0, 1]
R[my.data[,'IC'] == 0] <- my.data[my.data[,'IC'] == 0, 1]
L[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1, "L"]
R[my.data[,'IC'] == 1] <- my.data[my.data[,'IC'] == 1, "R"]
R[my.data[,"y"] == 2] <- Inf

TRT_C <- my.data[, "TRT_C"]
D <- cbind(L, R, TRT_C)

tm_wild <- system.time({
  fit1 <- ic_sp(cbind(D[,"L"], D[,"R"]) ~ TRT_C,
                model = 'po',
                data = as.data.frame(D),
                bs_samples = 100)
  s1 <- summary(fit1)
})
print(s1)
print_timing("KRAS wild-type (KRAS_C==0)", tm_wild)

#################### Total timing ####################
t_total <- proc.time() - t_total_start
print_timing("Total pipeline", t_total)
