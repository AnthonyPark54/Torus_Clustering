library(mvtnorm)
library(polynom)
library(igraph)

##
splitdata <- function(x){
  # Function for spliting the data
  n <- nrow(x)
  index <- sample(n, floor(n/2))
  x1 <- x[index,]
  x2 <- x[-index,]
  
  return(list(mix = x1, level = x2))
}

##
cal.parameter <- function(x1, ncomp = 3){
  # Function for obtaining parameter estimates for the bivariate von Mises Mixture
  summ <- summary(fit_vmsinmix(x1, method = "rwmh", ncomp = ncomp, unimodal.component = TRUE, n.iter = 1000))
  return(summ$`estimate `$mode)
}

##
conf.score.bvn.approx <- function(x, parameter){
  # Function for calculation of conformity scores using maximum of bivariate normal approximation
  k <- ncol(parameter)
  n <- nrow(x)
  
  x.shift1 <- x
  x.shift2 <- x
  x.shift3 <- x
  p.hat <- matrix(0, nrow = nrow(x), ncol = 4*k)
  conf.score.matrix <- matrix(0, nrow = nrow(x), ncol = k)
  for (i in 1:n){
    if (x[i,1] >= pi & x[i,2] >= pi){
      x.shift1[i,1] <- x[i,1] - 2*pi
      x.shift2[i,2] <- x[i,2] - 2*pi
      x.shift3[i,1] <- x[i,1] - 2*pi
      x.shift3[i,2] <- x[i,2] - 2*pi
    }else if (x[i,1] < pi & x[i,2] >= pi){
      x.shift1[i,1] <- x[i,1] + 2*pi
      x.shift2[i,2] <- x[i,2] - 2*pi
      x.shift3[i,1] <- x[i,1] + 2*pi
      x.shift3[i,2] <- x[i,2] - 2*pi
    }else if (x[i,1] < pi & x[i,2] < pi){
      x.shift1[i,1] <- x[i,1] + 2*pi
      x.shift2[i,2] <- x[i,2] + 2*pi
      x.shift3[i,1] <- x[i,1] + 2*pi
      x.shift3[i,2] <- x[i,2] + 2*pi
    }else{
      x.shift1[i,1] <- x[i,1] - 2*pi
      x.shift2[i,2] <- x[i,2] + 2*pi
      x.shift3[i,1] <- x[i,1] - 2*pi
      x.shift3[i,2] <- x[i,2] + 2*pi
    }
  }
  for (g in 1:k){
    mu1 <- parameter[5,g]
    mu2 <- parameter[6,g]
    sigma1 <- sqrt(parameter[3,g] / (parameter[2,g]*parameter[3,g] - parameter[4,g]^2))
    sigma2 <- sqrt(parameter[2,g] / (parameter[2,g]*parameter[3,g] - parameter[4,g]^2))
    rho <- parameter[4,g] / sqrt(parameter[2,g] * parameter[3,g])
    Sigma.mat <- matrix(c(sigma1^2, rho * sigma1 * sigma2, rho * sigma1 * sigma2, sigma2^2), nrow = 2, ncol = 2)
    
    p.hat[,g] <- log(parameter[1,g]) + log(dmvnorm(x, mean = c(mu1, mu2), sigma = Sigma.mat))
    p.hat[,g+k] <- log(parameter[1,g]) + log(dmvnorm(x.shift1, mean = c(mu1, mu2), sigma = Sigma.mat))
    p.hat[,g+2*k] <- log(parameter[1,g]) + log(dmvnorm(x.shift2, mean = c(mu1, mu2), sigma = Sigma.mat))
    p.hat[,g+3*k] <- log(parameter[1,g]) + log(dmvnorm(x.shift3, mean = c(mu1, mu2), sigma = Sigma.mat))
    conf.score.matrix[,g] <- apply(p.hat[,c(g,g+k,g+2*k,g+3*k)], 1, max)
  }
  
  conf.score <- apply(conf.score.matrix, 1, max)
  
  return(list(conf.score.matrix = conf.score.matrix, conf.score = conf.score))
}

##
conf.pred.set <- function(x, conf.score, parameter, alpha, grid.size = 200){
  n <- length(conf.score)
  ialpha <- floor((n+1)*alpha)
  level <- conf.score[order(conf.score)][ialpha]
  level.point <- x[order(conf.score),][ialpha,]
  
  Axis <- seq(0, 2*pi, length = grid.size)
  lattice <- cbind(rep(Axis, grid.size), rep(Axis, each = grid.size))
  
  y <- conf.score.bvn.approx(lattice, parameter)$conf.score
  conf.pred.set <- lattice[which(y >= level),]
  
  conf.pred.set[conf.pred.set[,1] < 0,] <- conf.pred.set[conf.pred.set[,1] < 0,] + 2*pi
  conf.pred.set[conf.pred.set[,1] > 2*pi,] <- conf.pred.set[conf.pred.set[,1] > 2*pi,] - 2*pi
  conf.pred.set[conf.pred.set[,2] < 0,] <- conf.pred.set[conf.pred.set[,2] < 0,] + 2*pi
  conf.pred.set[conf.pred.set[,2] > 2*pi,] <- conf.pred.set[conf.pred.set[,2] > 2*pi,] - 2*pi  
  
  return(list(Point = level.point, Level = level, Set = conf.pred.set))
}

##
Test.intersection.ellipse <- function(param.ellipse.1, param.ellipse.2, conf.score.level){
  prob.mix.1 <- param.ellipse.1[1]
  prob.mix.2 <- param.ellipse.2[1]
  
  mean.1 <- matrix(param.ellipse.1[5:6], ncol = 1)
  mean.2 <- matrix(param.ellipse.2[5:6], ncol = 1)
  
  sigma1.1 <- sqrt(param.ellipse.1[3] / (param.ellipse.1[2]*param.ellipse.1[3] - param.ellipse.1[4]^2))
  sigma2.1 <- sqrt(param.ellipse.1[2] / (param.ellipse.1[2]*param.ellipse.1[3] - param.ellipse.1[4]^2))
  rho.1 <- param.ellipse.1[4] / sqrt(param.ellipse.1[2] * param.ellipse.1[3])
  Sigma.mat.1 <- matrix(c(sigma1.1^2, rho.1 * sigma1.1 * sigma2.1, rho.1 * sigma1.1 * sigma2.1, sigma2.1^2), nrow = 2, ncol = 2)
  
  sigma1.2 <- sqrt(param.ellipse.2[3] / (param.ellipse.2[2]*param.ellipse.2[3] - param.ellipse.2[4]^2))
  sigma2.2 <- sqrt(param.ellipse.2[2] / (param.ellipse.2[2]*param.ellipse.2[3] - param.ellipse.2[4]^2))
  rho.2 <- param.ellipse.2[4] / sqrt(param.ellipse.2[2] * param.ellipse.2[3])
  Sigma.mat.2 <- matrix(c(sigma1.2^2, rho.2 * sigma1.2 * sigma2.2, rho.2 * sigma1.2 * sigma2.2, sigma2.2^2), nrow = 2, ncol = 2)
  
  Const.1 <- log(prob.mix.1) - log(2*pi*sigma1.1*sigma2.1*sqrt(1-rho.1^2))
  Const.2 <- log(prob.mix.2) - log(2*pi*sigma1.2*sigma2.2*sqrt(1-rho.2^2))
  
  M.1 <- solve(Sigma.mat.1) / (-2 * (conf.score.level - Const.1))
  M.2 <- solve(Sigma.mat.2) / (-2 * (conf.score.level - Const.2))
  
  # Now, each ellipse satisfies (x - mean.j)^T M.j (x - mean.j) = 1
  
  Eig.1 <- eigen(M.1)
  R.1 <- Eig.1$vectors
  D.1 <- diag((Eig.1$values))
  Eig.2 <- eigen(M.2)
  R.2 <- Eig.2$vectors
  D.2 <- diag((Eig.2$values))
  
  # Reduction to the unit circle and axis-aligned ellipse
  K.3 <- sqrt(D.1) %*% t(R.1) %*% (mean.2 - mean.1)
  M.3 <- diag((1/sqrt(diag(D.1)))) %*% t(R.1) %*% R.2 %*% D.2 %*% t(R.2) %*% R.1 %*% diag((1/sqrt(diag(D.1))))
  
  Eig.3 <- eigen(M.3)
  R <- Eig.3$vectors[,order(Eig.3$values)]
  D <- diag(Eig.3$values[order(Eig.3$values)])
  K <- t(R) %*% K.3
  
  # package 'polynom' is required
  d0 <- D[1,1]
  d1 <- D[2,2]
  k0 <- K[1]
  k1 <- K[2]
  coef.0 <- 1 - d0 * k0^2 - d1 * k1^2
  coef.1 <- -2 * d0 - 2 * d1 + 2 * d0 * d1 * k0^2 + 2 * d0 * d1 * k1^2
  coef.2 <- d0^2 + d1^2 + 4*d0*d1 - d0*d1^2*k0^2 - d0^2*d1*k1^2
  coef.3 <- -2 * d0 * d1^2 - 2 * d0^2 * d1
  coef.4 <- d0^2 * d1^2
  f.s <- polynomial(c(coef.0, coef.1, coef.2, coef.3, coef.4))
  root.f.s <- Re(solve(f.s))
  
  P.0 <- c(d0 * k0 * min(root.f.s) / (d0 * min(root.f.s) - 1), d1 * k1 * min(root.f.s) / (d1 * min(root.f.s) - 1))
  P.1 <- c(d0 * k0 * max(root.f.s) / (d0 * max(root.f.s) - 1), d1 * k1 * max(root.f.s) / (d1 * max(root.f.s) - 1))
  
  # From now on, test whether two ellipses are overlapped or not
  minDistance <- sqrt(sum(P.0^2))
  maxDistance <- sqrt(sum(P.1^2))
  Ind.Overlap <- 0
  if (maxDistance <= 1){
    Ind.Overlap <- 1
  }else{
    if (minDistance < 1){
      Ind.Overlap <- 1
    }else if (minDistance > 1){
      if (d0 * k0^2 + d1 * k1^2 - 1 > 0){
        Ind.Overlap <- 0
      }else{
        Ind.Overlap <- 1
      }
    }else{
      if (d0 * k0^2 + d1 * k1^2 - 1 > 0){
        Ind.Overlap <- 0
      }else{
        Ind.Overlap <- 1
      }
    }
  }
  return(Ind.Overlap)
}
conf.pred.clusters <- function(parameter, conf.score.level){
  K <- ncol(parameter)
  
  Matrix.Cluster.Indicator <- matrix(0, ncol = K, nrow = K)
  for (rr in 1:(K-1)){
    for (cc in (rr+1):K){
      Test.rr.cc.temp <- Test.intersection.ellipse(parameter[,rr], parameter[,cc], conf.score.level = conf.score.level)
      if (Test.rr.cc.temp == 0){
        parameter.shift <- matrix(rep(parameter[,cc],3), ncol = 3, byrow = FALSE)
        if (parameter[5,cc] >= pi & parameter[6,cc] >= pi){
          parameter.shift[5,1] <- parameter.shift[5,1] - 2*pi
          parameter.shift[6,2] <- parameter.shift[6,2] - 2*pi
          parameter.shift[5,3] <- parameter.shift[5,3] - 2*pi
          parameter.shift[6,3] <- parameter.shift[6,3] - 2*pi
        }else if (parameter[5,cc] < pi & parameter[6,cc] >= pi){
          parameter.shift[5,1] <- parameter.shift[5,1] + 2*pi
          parameter.shift[6,2] <- parameter.shift[6,2] - 2*pi
          parameter.shift[5,3] <- parameter.shift[5,3] + 2*pi
          parameter.shift[6,3] <- parameter.shift[6,3] - 2*pi
        }else if (parameter[5,cc] < pi & parameter[6,cc] < pi){
          parameter.shift[5,1] <- parameter.shift[5,1] + 2*pi
          parameter.shift[6,2] <- parameter.shift[6,2] + 2*pi
          parameter.shift[5,3] <- parameter.shift[5,3] + 2*pi
          parameter.shift[6,3] <- parameter.shift[6,3] + 2*pi
        }else{
          parameter.shift[5,1] <- parameter.shift[5,1] - 2*pi
          parameter.shift[6,2] <- parameter.shift[6,2] + 2*pi
          parameter.shift[5,3] <- parameter.shift[5,3] - 2*pi
          parameter.shift[6,3] <- parameter.shift[6,3] + 2*pi
        }
        Test.rr.cc.temp <- Test.intersection.ellipse(parameter[,rr], parameter.shift[,1], conf.score.level = conf.score.level)
        if (Test.rr.cc.temp == 0){
          Test.rr.cc.temp <- Test.intersection.ellipse(parameter[,rr], parameter.shift[,2], conf.score.level = conf.score.level)
        }
        if (Test.rr.cc.temp == 0){
          Test.rr.cc.temp <- Test.intersection.ellipse(parameter[,rr], parameter.shift[,3], conf.score.level = conf.score.level)
        }
      }
      Matrix.Cluster.Indicator[rr,cc] <- Test.rr.cc.temp
      Matrix.Cluster.Indicator[cc,rr] <- Matrix.Cluster.Indicator[rr,cc]
    }
  }
  Graph <- graph_from_adjacency_matrix(Matrix.Cluster.Indicator)
  Components <- components(Graph, mode = "weak")$membership
  Index.Clusters <- data.frame(K = 1:K,
                               Cluster = Components)
  return(Index.Clusters)
}

##
conf.pred.assign.clusters <- function(x, parameter, conf.score.level, Index.Clusters, option.clusters = 1){
  if (option.clusters == 1){
    x <- as.matrix(x)
    K <- ncol(parameter)
    N <- nrow(x)
    Distance.Mat <- matrix(1e5, ncol = K, nrow = N)
    Cluster.Assignment <- data.frame(x = x,
                                     Cluster = rep(1, N))
    for (cc in 1:K){
      Mu <- matrix(parameter[5:6,cc], ncol = 1)
      sigma1 <- sqrt(parameter[3,cc] / (parameter[2,cc]*parameter[3,cc] - parameter[4,cc]^2))
      sigma2 <- sqrt(parameter[2,cc] / (parameter[2,cc]*parameter[3,cc] - parameter[4,cc]^2))
      rho <- parameter[4,cc] / sqrt(parameter[2,cc] * parameter[3,cc])
      Sigma.mat <- matrix(c(sigma1^2, rho * sigma1 * sigma2, rho * sigma1 * sigma2, sigma2^2), nrow = 2, ncol = 2)
      for (rr in 1:N){
        Distance.Mat[rr,cc] <- sqrt((x[rr,] - t(Mu)) %*% solve(Sigma.mat) %*% t(x[rr,] - t(Mu)))
        if (x[rr,1] >= pi & x[rr,2] >= pi){
          Distance.temp <- min(sqrt(((x[rr,] + matrix(c(-2*pi,0),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(-2*pi,0),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(0,-2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(0,-2*pi),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(-2*pi,-2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(-2*pi,-2*pi),nrow=1)) - t(Mu))))
        }else if (x[rr,1] < pi & x[rr,2] >= pi){
          Distance.temp <- min(sqrt(((x[rr,] + matrix(c(2*pi,0),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(2*pi,0),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(0,-2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(0,-2*pi),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(2*pi,-2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(2*pi,-2*pi),nrow=1)) - t(Mu))))
        }else if (x[rr,1] < pi & x[rr,2] < pi){
          Distance.temp <- min(sqrt(((x[rr,] + matrix(c(2*pi,0),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(2*pi,0),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(0,2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(0,2*pi),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(2*pi,2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(2*pi,2*pi),nrow=1)) - t(Mu))))
        }else{
          Distance.temp <- min(sqrt(((x[rr,] + matrix(c(-2*pi,0),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(-2*pi,0),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(0,2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(0,2*pi),nrow=1)) - t(Mu))),
                               sqrt(((x[rr,] + matrix(c(-2*pi,2*pi),nrow=1)) - t(Mu)) %*% solve(Sigma.mat) %*% t((x[rr,] + matrix(c(-2*pi,2*pi),nrow=1)) - t(Mu))))
        }
        Distance.Mat[rr,cc] <- min(Distance.Mat[rr,cc], Distance.temp)
        Cluster.Assignment[rr,3] <- Index.Clusters$Cluster[which.min(Distance.Mat[rr,])]
      }
    }
    return(list(Mahal.dist = Distance.Mat, Cluster.Assignment = Cluster.Assignment))
  }else if (option.clusters == 2){
    N <- nrow(x)
    conf.score.matrix <- conf.score.bvn.approx(x, parameter)$conf.score.matrix
    
    Cluster.Assignment <- data.frame(x = x, Cluster = rep(0, N))
    Cluster.Assignment$Cluster[apply(conf.score.matrix, 1, max) >= conf.score.level] <- Index.Clusters$Cluster[apply(conf.score.matrix[apply(conf.score.matrix, 1, max) >= conf.score.level,], 1, which.max)]
    return(list(Conf.Scores = conf.score.matrix, Cluster.Assignment = Cluster.Assignment))
  }else{
    print("Not available option!")
  }
}

##
Draw.ellipses.bvn.approx <- function(data, parameter.ellipse, conf.score.level, cut.point, assigned.clusters){
  par(pty = "s")
  plot(0,0, type = "n", xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "", main ="", xaxt = "n", yaxt = "n")
  k <- ncol(parameter.ellipse)
  for (g in 1:k){
    mu1 <- parameter.ellipse[5,g]
    mu2 <- parameter.ellipse[6,g]
    sigma1 <- sqrt(parameter.ellipse[3,g] / (parameter.ellipse[2,g]*parameter.ellipse[3,g] - parameter.ellipse[4,g]^2))
    sigma2 <- sqrt(parameter.ellipse[2,g] / (parameter.ellipse[2,g]*parameter.ellipse[3,g] - parameter.ellipse[4,g]^2))
    rho <- parameter.ellipse[4,g] / sqrt(parameter.ellipse[2,g] * parameter.ellipse[3,g])
    Sigma.mat <- matrix(c(sigma1^2, rho * sigma1 * sigma2, rho * sigma1 * sigma2, sigma2^2), nrow = 2, ncol = 2)
    
    Const <- log(parameter.ellipse[1,g]) - log(2*pi*sigma1*sigma2*sqrt(1-rho^2))
    
    Sigma.inv.tilde <- solve(Sigma.mat) / (-2 * (conf.score.level - Const))
    
    Eig.Sigma.inv.tilde <- eigen(Sigma.inv.tilde)
    R <- Eig.Sigma.inv.tilde$vectors
    D <- diag(Eig.Sigma.inv.tilde$values)
    
    theta <- seq(0, 2*pi, length.out = 1000)
    x <- cos(theta)
    y <- sin(theta)
    Z <- cbind(x,y)
    Y1 <- Z
    Y2 <- Z
    Y3 <- Z
    Y4 <- Z
    for (i in 1:1000){
      Y1[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1, mu2)
      if (mu1 >= pi & mu2 >= pi){
        Y2[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 - 2*pi, mu2)
        Y3[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1, mu2 - 2*pi)
        Y4[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 - 2*pi, mu2 - 2*pi)
      }else if (mu1 < pi & mu2 >= pi){
        Y2[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2)
        Y3[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1, mu2 - 2*pi)
        Y4[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2 - 2*pi)
      }else if (mu1 < pi & mu2 < pi){
        Y2[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2)
        Y3[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1, mu2 + 2*pi)
        Y4[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2 + 2*pi)
      }else{
        Y2[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2)
        Y3[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1, mu2 - 2*pi)
        Y4[i,] <- R %*% diag(1/sqrt(diag(D))) %*% Z[i,] + c(mu1 + 2*pi, mu2 - 2*pi)
      }
    }
    polygon(Y1[,1], Y1[,2], xlim = c(0,2*pi), ylim = c(0,2*pi),
            border = NA, col = rgb(col2rgb("grey")[1], col2rgb("grey")[2], col2rgb("grey")[3], max = 255, alpha = 80 * 255 / 100))
    polygon(Y2[,1], Y2[,2], xlim = c(0,2*pi), ylim = c(0,2*pi),
            border = NA, col = rgb(col2rgb("grey")[1], col2rgb("grey")[2], col2rgb("grey")[3], max = 255, alpha = 80 * 255 / 100))
    polygon(Y3[,1], Y3[,2], xlim = c(0,2*pi), ylim = c(0,2*pi),
            border = NA, col = rgb(col2rgb("grey")[1], col2rgb("grey")[2], col2rgb("grey")[3], max = 255, alpha = 80 * 255 / 100))
    polygon(Y4[,1], Y4[,2], xlim = c(0,2*pi), ylim = c(0,2*pi),
            border = NA, col = rgb(col2rgb("grey")[1], col2rgb("grey")[2], col2rgb("grey")[3], max = 255, alpha = 80 * 255 / 100))
    lines(Y1[,1], Y1[,2], type = "l", xlim = c(0,2*pi), ylim = c(0,2*pi), col = "black")
    lines(Y2[,1], Y2[,2], type = "l", xlim = c(0,2*pi), ylim = c(0,2*pi), col = "black")
    lines(Y3[,1], Y3[,2], type = "l", xlim = c(0,2*pi), ylim = c(0,2*pi), col = "black")
    lines(Y4[,1], Y4[,2], type = "l", xlim = c(0,2*pi), ylim = c(0,2*pi), col = "black")
    points(mu1, mu2,
           type = "p", pch = 15, cex = 1.5, col = "red",
           xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "")
  }
  text(data[,1], data[,2], label = as.character(assigned.clusters), cex = 0.7, xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "")
  points(cut.point[1], cut.point[2],
         type = "p", pch = "*", cex = 3, col = "blue", 
         xlim = c(0,2*pi), ylim = c(0,2*pi), xlab = "", ylab = "")
  axis(side = 1, at = c(0, pi/2, pi, 3/2*pi, 2*pi),
       labels = c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)),
       las = 1, cex.axis = 1.1)
  axis(side = 2, at = c(0, pi/2, pi, 3/2*pi, 2*pi),
       labels = c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)),
       las = 2, cex.axis = 1.1)
}
