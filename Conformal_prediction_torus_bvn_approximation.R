library(bio3d)
library(MASS)
library(gsl)
library(BAMBI)
library(mvtnorm)
library(polynom)
library(igraph)

source("Conformal_Prediction_Torus_BVNapprox_functions.R")

set.seed(12345)
## Trivial example without overlapping
G1 <- mvrnorm(n = 100, mu = c(pi/2,pi), Sigma = matrix(c(0.25,0,0,0.25),nrow = 2))
G2 <- mvrnorm(n = 100, mu = c(3*pi/2,pi), Sigma = matrix(c(0.25,0,0,0.25),nrow = 2))
plot(G1, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue")
lines(G2, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "red")

Example1 <- rbind(G1, G2)

## Trivial example with overlapping
G1 <- mvrnorm(n = 100, mu = c(3*pi/4,pi), Sigma = matrix(c(0.3,0,0,0.3),nrow = 2))
G2 <- mvrnorm(n = 100, mu = c(5*pi/4,pi), Sigma = matrix(c(0.3,0,0,0.3),nrow = 2))
plot(G1, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "blue")
lines(G2, type = "p", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", col = "red")

Example2 <- rbind(G1, G2)

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

# Step 1: split the data
Example1.split <- splitdata(Example1) 
Example2.split <- splitdata(Example2)
Example3.split <- splitdata(Example3)
Example4.split <- splitdata(Example4)

# Step 2: Fit the bivariate von Mises (sine model) mixture with known K = 2
fit.vMM.parameter.Ex1 <- cal.parameter(Example1.split$mix, ncomp = 2)
fit.vMM.parameter.Ex2 <- cal.parameter(Example2.split$mix, ncomp = 2)
fit.vMM.parameter.Ex3 <- cal.parameter(Example3.split$mix, ncomp = 2)
fit.vMM.parameter.Ex4 <- cal.parameter(Example4.split$mix, ncomp = 7)

# Step 3: Set the conformity score as max of approximated bivariate normal
Conf.Scores.Ex1 <- conf.score.bvn.approx(Example1.split$level, fit.vMM.parameter.Ex1)
Conf.Scores.Ex2 <- conf.score.bvn.approx(Example2.split$level, fit.vMM.parameter.Ex2)
Conf.Scores.Ex3 <- conf.score.bvn.approx(Example3.split$level, fit.vMM.parameter.Ex3)
Conf.Scores.Ex4 <- conf.score.bvn.approx(Example4.split$level, fit.vMM.parameter.Ex4)

# Step 4: Obtain conformal prediction set for given alpha
Conf.Pred.Set.Ex1 <- conf.pred.set(Example1.split$level, Conf.Scores.Ex1$conf.score, 
                                   fit.vMM.parameter.Ex1, alpha = 0.05, grid.size = 200)
Conf.Pred.Set.Ex2 <- conf.pred.set(Example2.split$level, Conf.Scores.Ex2$conf.score, 
                                   fit.vMM.parameter.Ex2, alpha = 0.05, grid.size = 200)
Conf.Pred.Set.Ex3 <- conf.pred.set(Example3.split$level, Conf.Scores.Ex3$conf.score, 
                                   fit.vMM.parameter.Ex3, alpha = 0.05, grid.size = 200)
Conf.Pred.Set.Ex4 <- conf.pred.set(Example4.split$level, Conf.Scores.Ex4$conf.score,
                                   fit.vMM.parameter.Ex4, alpha = 0.05, grid.size = 200)

# Step 5: Determine clusters from the ellipses
Index.Clusters.Ex1 <- conf.pred.clusters(parameter = fit.vMM.parameter.Ex1, conf.score.level = Conf.Pred.Set.Ex1$Level)
Index.Clusters.Ex2 <- conf.pred.clusters(parameter = fit.vMM.parameter.Ex2, conf.score.level = Conf.Pred.Set.Ex2$Level)
Index.Clusters.Ex3 <- conf.pred.clusters(parameter = fit.vMM.parameter.Ex3, conf.score.level = Conf.Pred.Set.Ex3$Level)
Index.Clusters.Ex4 <- conf.pred.clusters(parameter = fit.vMM.parameter.Ex4, conf.score.level = Conf.Pred.Set.Ex4$Level)

Example1.Cluster.Assignment <- conf.pred.assign.clusters(x = Example1, parameter = fit.vMM.parameter.Ex1,
                                                         conf.score.level = Conf.Pred.Set.Ex1$Level,
                                                         Index.Clusters = Index.Clusters.Ex1,
                                                         option.clusters = 1)
Example2.Cluster.Assignment <- conf.pred.assign.clusters(x = Example2, parameter = fit.vMM.parameter.Ex2,
                                                         conf.score.level = Conf.Pred.Set.Ex2$Level,
                                                         Index.Clusters = Index.Clusters.Ex2,
                                                         option.clusters = 1)
Example3.Cluster.Assignment <- conf.pred.assign.clusters(x = Example3, parameter = fit.vMM.parameter.Ex3,
                                                         conf.score.level = Conf.Pred.Set.Ex3$Level,
                                                         Index.Clusters = Index.Clusters.Ex3,
                                                         option.clusters = 1)
Example4.Cluster.Assignment <- conf.pred.assign.clusters(x = Example4, parameter = fit.vMM.parameter.Ex4,
                                                         conf.score.level = Conf.Pred.Set.Ex4$Level,
                                                         Index.Clusters = Index.Clusters.Ex4,
                                                         option.clusters = 1)
#--------------------------------------------------------------------------------------------------------------
Example1.Cluster.Assignment.op2 <- conf.pred.assign.clusters(x = Example1, parameter = fit.vMM.parameter.Ex1,
                                                         conf.score.level = Conf.Pred.Set.Ex1$Level,
                                                         Index.Clusters = Index.Clusters.Ex1,
                                                         option.clusters = 2)
Example2.Cluster.Assignment.op2 <- conf.pred.assign.clusters(x = Example2, parameter = fit.vMM.parameter.Ex2,
                                                         conf.score.level = Conf.Pred.Set.Ex2$Level,
                                                         Index.Clusters = Index.Clusters.Ex2,
                                                         option.clusters = 2)
Example3.Cluster.Assignment.op2 <- conf.pred.assign.clusters(x = Example3, parameter = fit.vMM.parameter.Ex3,
                                                         conf.score.level = Conf.Pred.Set.Ex3$Level,
                                                         Index.Clusters = Index.Clusters.Ex3,
                                                         option.clusters = 2)
Example4.Cluster.Assignment.op2 <- conf.pred.assign.clusters(x = Example4, parameter = fit.vMM.parameter.Ex4,
                                                         conf.score.level = Conf.Pred.Set.Ex4$Level,
                                                         Index.Clusters = Index.Clusters.Ex4,
                                                         option.clusters = 2)
#--------------------------------------------------------------------------------------------------------------

# Step 6 (optional): Draw fitted ellipses with conformal prediction boundaries
Draw.ellipses.bvn.approx(data = Example1, parameter.ellipse = fit.vMM.parameter.Ex1, 
                         conf.score.level = Conf.Pred.Set.Ex1$Level,
                         cut.point = Conf.Pred.Set.Ex1$Point,
                         assigned.clusters = Example1.Cluster.Assignment$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example2, parameter.ellipse = fit.vMM.parameter.Ex2, 
                         conf.score.level = Conf.Pred.Set.Ex2$Level,
                         cut.point = Conf.Pred.Set.Ex2$Point,
                         assigned.clusters = Example2.Cluster.Assignment$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example3, parameter.ellipse = fit.vMM.parameter.Ex3, 
                         conf.score.level = Conf.Pred.Set.Ex3$Level,
                         cut.point = Conf.Pred.Set.Ex3$Point,
                         assigned.clusters = Example3.Cluster.Assignment$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example4, parameter.ellipse = fit.vMM.parameter.Ex4, 
                         conf.score.level = Conf.Pred.Set.Ex4$Level,
                         cut.point = Conf.Pred.Set.Ex4$Point,
                         assigned.clusters = Example4.Cluster.Assignment$Cluster.Assignment[,3])

#--------------------------------------------------------------------------------------------------------------
Draw.ellipses.bvn.approx(data = Example1, parameter.ellipse = fit.vMM.parameter.Ex1, 
                         conf.score.level = Conf.Pred.Set.Ex1$Level,
                         cut.point = Conf.Pred.Set.Ex1$Point,
                         assigned.clusters = Example1.Cluster.Assignment.op2$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example2, parameter.ellipse = fit.vMM.parameter.Ex2, 
                         conf.score.level = Conf.Pred.Set.Ex2$Level,
                         cut.point = Conf.Pred.Set.Ex2$Point,
                         assigned.clusters = Example2.Cluster.Assignment.op2$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example3, parameter.ellipse = fit.vMM.parameter.Ex3, 
                         conf.score.level = Conf.Pred.Set.Ex3$Level,
                         cut.point = Conf.Pred.Set.Ex3$Point,
                         assigned.clusters = Example3.Cluster.Assignment.op2$Cluster.Assignment[,3])
Draw.ellipses.bvn.approx(data = Example4, parameter.ellipse = fit.vMM.parameter.Ex4, 
                         conf.score.level = Conf.Pred.Set.Ex4$Level,
                         cut.point = Conf.Pred.Set.Ex4$Point,
                         assigned.clusters = Example4.Cluster.Assignment.op2$Cluster.Assignment[,3])
#--------------------------------------------------------------------------------------------------------------





