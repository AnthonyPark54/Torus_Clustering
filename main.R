
# 2019 Research SNU
# R code by Kiho Park with Prof. Sung kyu Jung

##### Bio Data ##############################################
# https://proteinstructures.com/Structure/Structure/Ramachandran-plot.html
# choose in http://www.rcsb.org

library(MASS)
library(bio3d)
source("vMM.R")
source("kernel.R")

pdb <- read.pdb("6a32")
a <- torsion.pdb(pdb)
data <- cbind(a$phi/180*pi,a$psi/180*pi)
data <- data[-which(is.na(data[,1])|is.na(data[,2])),]

on.torus <- function(x){
  # range to [0,2*pi) x [0,2*pi)
  
  y <- x
  y[,1] <- x[,1] - 2*pi*floor(x[,1]/2/pi)
  y[,2] <- x[,2] - 2*pi*floor(x[,2]/2/pi)
  return(y)
}

data <- on.torus(data)

##### New Data ##############################################
# Example: ball and L shape
Mu <- c(1,5)
Sigma <- matrix(c(0.05,0,0,0.05),2,2)

unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)

Example1 <- rbind(mvrnorm(n=100, Mu, Sigma), unidata.unif, data.unif1, data.unif2)
Example1 <- on.torus(Example1)
plot(Example1, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xlab="phi", ylab="psi")

# Example: five clusters
Mu1 <- c(3,0)
Mu2 <- c(2,2)
Mu3 <- c(1,4)
Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)

unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4), sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))

Example2 <- rbind(mvrnorm(n=70, Mu1, Sigma1), mvrnorm(n=50, Mu2, Sigma2), mvrnorm(n=50, Mu3, Sigma3), data.unif, data.diamond)
Example2 <- on.torus(Example2)
plot(Example2, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xlab="phi", ylab="psi")

## Example: mean on the borderline
G1 <- mvrnorm(n = 100, mu = c(0,15*pi/8), Sigma = matrix(c(0.3,0,0,0.3),nrow = 2))
G1[G1[,1] < 0,1] <- G1[G1[,1] < 0,1] + 2*pi
G1[G1[,2] > 2*pi,2] <- G1[G1[,2] > 2*pi,2] - 2*pi
G2 <- mvrnorm(n = 100, mu = c(pi,15*pi/8), Sigma = matrix(c(0.3,0,0,0.3),nrow = 2))
G2[G2[,2] > 2*pi,2] <- G2[G2[,2] > 2*pi,2] - 2*pi
plot(G1, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue")
lines(G2, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "red")

Example3 <- rbind(G1, G2)

## Example: Complex (several clusters and some are on borderline)
G1 <- mvrnorm(n = 50, mu = c(1/10*2*pi, 1/10*2*pi), Sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))
plot(G1, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue", pch = 1)
G2 <- mvrnorm(n = 50, mu = c(3/10*2*pi, 1/10*2*pi), Sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))
points(G2, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue", pch = 2)
G3 <- mvrnorm(n = 50, mu = c(1/10*2*pi, 3/10*2*pi), Sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))
points(G3, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue", pch = 3)
G4 <- mvrnorm(n = 50, mu = c(5/10*2*pi, 1/10*2*pi), Sigma = matrix(c(0.05, 0, 0, 0.05), nrow = 2))
points(G4, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue", pch = 4)
G5 <- mvrnorm(n = 50, mu = c(1/10*2*pi, 5/10*2*pi), Sigma = matrix(c(0.05, 0, 0, 0.05), nrow = 2))
points(G5, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue", pch = 5)
G6 <- mvrnorm(n = 50, mu = c(9.5/10*2*pi, 8/10*2*pi), Sigma = matrix(c(0.4, 0, 0, 0.05), nrow = 2))
G6[G6[,1] > 2*pi,1] <- G6[G6[,1] > 2*pi,1] - 2*pi
G6[G6[,2] > 2*pi,2] <- G6[G6[,2] > 2*pi,2] - 2*pi
points(G6, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "red", pch = 6)
G7 <- mvrnorm(n = 50, mu = c(8/10*2*pi, 4/10*2*pi), Sigma = matrix(c(0.05, 0, 0, 0.05), nrow = 2))
points(G7, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "green", pch = 7)

Example4 <- rbind(G1, G2, G3, G4, G5, G6, G7)


##################### Application #############

### Kernel Density Method
Set.Kernel <- kernel.method(data, alpha = 0.1,  kap1 = 50, kap2 = 50, grid.size = 100, conformal.set = TRUE)
Set.Kernel$conformal.plot

# KDE optimization
original.ker.optimal <- kernel.optimization(data)
ex1.ker.optimal <- kernel.optimization(Example1)
ex2.ker.optimal <- kernel.optimization(Example2)
ex3.ker.optimal <- kernel.optimization(Example3)
ex4.ker.optimal <- kernel.optimization(Example4)

### von Mises mixture method
Parm.vMM <- vMM.method.parameter(data, ncomp = 5)
Set.vMM <- vMM.method.alpha(first.cal = Parm.vMM, alpha = 0.05, grid.size = 100)
Set.vMM$total.plot

# vMM optimization
original.vMM.optimal <- vMM.optimization(data)
ex1.vMM.optimal <- vMM.optimization(Example1)
ex2.vMM.optimal <- vMM.optimization(Example2)
ex3.vMM.optimal <- vMM.optimization(Example3)
ex4.vMM.optimal <- vMM.optimization(Example4)

Clustering.vMM(original.vMM.optimal)
Clustering.vMM(ex1.vMM.optimal)
Clustering.vMM(ex2.vMM.optimal)
Clustering.vMM(ex3.vMM.optimal)
Clustering.vMM(ex4.vMM.optimal)

