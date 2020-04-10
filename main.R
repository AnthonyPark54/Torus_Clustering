
# 2019 Research SNU
# R code by Kiho Park with Prof. Sung kyu Jung

##### Bio Data ##############################################
# https://proteinstructures.com/Structure/Structure/Ramachandran-plot.html
# choose in http://www.rcsb.org
library(bio3d)
pdb <- read.pdb("6a32")
a <- torsion.pdb(pdb)
data <- cbind(a$phi/180*pi,a$psi/180*pi)
data <- data[-which(is.na(data[,1])|is.na(data[,2])),]

# [-pi,pi) x [-pi,pi) to [0,2*pi) x [0,2*pi)
on.torus <- function(x){
  y <- x
  y[,1] <- x[,1] - 2*pi*floor(x[,1]/2/pi)
  y[,2] <- x[,2] - 2*pi*floor(x[,2]/2/pi)
  return(y)
}

data <- on.torus(data)

### Kernel Density Method

Set.Kernel <- kernel.method(data, conformal.set = TRUE, alpha = 0.1,  kap1 = 50, kap2 = 50, div.num = 100)
Set.Kernel$conformal.plot
Optimal.Set.Kernel <- kernel.optimization(data, first.alpha = 0.1)

### von Mises mixture method

Parm.vMM <- vMM.method.parameter(data, ncomp = 5)
Set.vMM <- vMM.method.alpha(first.cal = Parm.vMM, alpha = 0.05, div.num = 100)
Set.vMM$total.plot
Optimal.Set.vMM <- vMM.optimization(data, first.alpha = 0.1)


##### New Data ##############################################
# New data 1
library(MASS)
Mu1 <- c(3,0)
Mu2 <- c(2,2)
Mu3 <- c(1,4)
Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)

unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4), sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))

new.data1 <- rbind(mvrnorm(n=70, Mu1, Sigma1), mvrnorm(n=50, Mu2, Sigma2), mvrnorm(n=50, Mu3, Sigma3), data.unif, data.diamond)
new.data1 <- on.torus(new.data1)

# New data 2
library(MASS)
Mu <- c(1,5)
Sigma <- matrix(c(0.05,0,0,0.05),2,2)

unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)

new.data2 <- rbind(mvrnorm(n=100, Mu, Sigma), unidata.unif, data.unif1, data.unif2)
new.data2 <- on.torus(new.data2)

# New data 3
library(MASS)
new.data3 <- rbind(mvrnorm(n=200, Mu, Sigma), data.unif1, data.unif2)
new.data3 <- on.torus(new.data3)

