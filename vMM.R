
# This function conducts conformal prediction sets on the torus data with von Mises mixture model.
# Finally, conduct C_alpha on the torus data.
# For clustering, this function find optimal prediction set.

###########################################################
## Torus Data : n x 2 dataframe ###########################
#### each row : [0,2*pi) x [0,2*pi) #######################
###########################################################
## Parameter setting ######################################
# alpha = 0.05
# ncomp = 25
# grid.size = 100
# first.alpha = 0.1
###########################################################

source("Conformal_Prediction_Torus_BVNapprox_functions.R")

### basic functions ###

splitdata <- function(data){
  # Split a data to apply the inductive conformal prediction
  # input : (n x 2 matrix) torus data
  # output : (list with elements mix two n/2 x 2 matrix) two randomly splited data sets
  
  n <- nrow(data)
  index <- sample(n, floor(n/2))
  x1 <- data[index,]
  x2 <- data[-index,]
  
  return(list(mix = x1, level = x2))
}

cal.parameter <- function(data, ncomp = 3){
  # Fit von Mises mixture model to a data
  # input : (n x 2 matrix) torus data, number of component K
  # output : (6 x K matrix) parameters of von Mises distribution
  
  library(BAMBI)
  summ <- summary(fit_vmsinmix(data, method = "rwmh", ncomp = ncomp, unimodal.component = TRUE, n.iter = 3000))
  return(summ$`estimate `$mode)
}

cal.constant <- function(parameter){
  # Calculate constants for the von Mises pdf 
  # input : (6 x K matrix) vMM parameter
  # output : (K vector) Constant for each component
  
  k <- ncol(parameter)
  d <- parameter
  
  grid.size <- 200
  Axis <- seq(0, 2*pi, length = grid.size)
  lattice <- cbind(rep(Axis, grid.size),rep(Axis, each=grid.size))
  
  const <- rep(0, k)
  
  y <- matrix(0, nrow(lattice), k)
  for(i in 1:k){
    for(j in 1:nrow(lattice)){
      y[j,i] <- exp(d[2,i]*cos(lattice[j,1]-d[5,i])
                    + d[3,i]*cos(lattice[j,2]-d[6,i])
                    + d[4,i]*sin(lattice[j,1]-d[5,i])*sin(lattice[j,2]-d[6,i]))
    }
    const[i] <- grid.size^2/sum(y[,i])/4/pi/pi
  }
  return(const)
}

cal.pdf.each <- function(lattice, parameter, constant){
  # Calculate the pdf of a point by each component
  # input : (m x 2 matrix) lattice, (6 x K matrix) vMM parameter, (K vector) constant for each component
  # output : (m x K matrix) each weighted pdf of a point
  
  k <- ncol(parameter)
  d <- parameter
  
  y <- matrix(0, nrow(lattice), k)
  
  for(i in 1:k){
    for(j in 1:nrow(lattice)){
      y[j,i] <- exp(d[2,i]*cos(lattice[j,1]-d[5,i])
                    + d[3,i]*cos(lattice[j,2]-d[6,i])
                    + d[4,i]*sin(lattice[j,1]-d[5,i])*sin(lattice[j,2]-d[6,i]))
    }
    y[,i] <- d[1,i]*y[,i]*constant[i]
  }
  return(y)
}

cal.level <- function(data, parameter, constant, alpha = 0.05, maximum = FALSE){
  # Calculate pdf level for alpha
  # input : (n x 2 matrix) data, (6 x K matrix) vMM parameter, (K vector) constant for each component, level alpha, level with maximum "TRUE"
  # output : pdf level
  
  pdf.each <- cal.pdf.each(data, parameter, constant)
  if(maximum){
    pdf <- apply(pdf.each, 1, max)
  } else{
    pdf <- rowSums(pdf.each)
  }
  
  n <- nrow(data)
  ialpha <- floor((n+1)*alpha)
  return(pdf[order(pdf)][ialpha])
}

# VMM
cal.conformal.set <- function(parameter, constant, pdf.level, grid.size = 100){
  # Obtain conformal set from vMM and pdf.level
  # input : (6 x K matrix) vMM parameter, (K vector) constant for each component, pdf level, grid size
  # output : (m x 3 matrix) conformal prediction set with pdf
  
  k <- ncol(parameter)
  d <- parameter
  names <- c("phi","psi", "pdf")
  
  Axis <- seq(0, 2*pi, length = grid.size)
  lattice <- cbind(rep(Axis, grid.size),rep(Axis, each=grid.size))
  
  y <- cal.pdf.each(lattice, parameter, constant)
  pdf <- rowSums(y)
  conf <- cbind(lattice[which(pdf >= pdf.level),], pdf[which(pdf >= pdf.level)])
  conf <- as.data.frame(conf)
  colnames(conf) <- names
  return(conf)
}

visual.vMM <- function(data, conf, ncomp = 3, alpha = 0.1){
  # Visualize conformal prediction set
  # input : (n x 2 matrix) data, (m x 3 matrix) conformal prediction set , number of component K, level alpha
  # output : conformal prediction set plot
  
  library(ggplot2)
  data_plot <- as.data.frame(data)
  colnames(data_plot) <- c("phi", "psi")
  
  gg <- ggplot(conf, aes(x=phi, y=psi)) +  geom_raster(aes(fill=pdf), interpolate = TRUE) +
    scale_fill_gradientn(colours=c("green","white"), guide = FALSE) +
    coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
    geom_point(aes(x=phi, y=psi), data_plot)  +
    # ggtitle(paste("Plot of conformal prediction set with vMM (K=",  ncomp, ", alpha=", alpha, ")", sep="")) +
    # theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(~ phi)) + ylab(expression(~ psi))
  return(gg)
}

# ellipse
cal.ellipse <- function(parameter, constant, pdf.level){
  # Approximate parameters for ellipse
  # input : (6 x K matrix) vMM parameter, (K vector) constant for each component, pdf level
  # output : (list) parameters of ellipses
  
  d <- parameter
  k <- ncol(d)
  
  mu.array <- matrix(0, nrow = 2, ncol = k)
  cov.array <- list()
  r.array <- c()
  for(i in 1:k){
    mu.array[,i] <- d[5:6, i]
    cov.array[[i]] <- solve(matrix(c(d[2,i]/2, -d[4,i]/2, -d[4,i]/2, d[3,i]/2), 2))
    r.array[i] <- sqrt(d[2,i]+d[3,i]-log(pdf.level)+log(constant[i])+log(d[1,i]))
  }
  return(list(mu = mu.array, cov.shape = cov.array, r = r.array))
}

visual.ellipse <- function(data, ellipse){
  # Visualize ellipses
  # input : (n x 2 matrix) data, (list) ellipse
  # output : ellipse plot
  
  library(ggplot2)
  library(ggforce)
  data_plot <- as.data.frame(data)
  colnames(data_plot) <- c("phi", "psi")
  
  k <- length(ellipse$r)
  
  color.array <- c("yellow", "pink", "green", "cyan", "salmon", "brwon", "cornsilk", "olivegreen", "chocolate", "purple")
  mu.x <- matrix()
  mu.y <- c()
  a <- c()
  b <- c()
  angle <- c()
  
  for(i in 1:k){
    mu.x[(9*i-8):(9*i-6)] <- ellipse$mu[1,i] - 2*pi
    mu.x[(9*i-5):(9*i-3)] <- ellipse$mu[1,i]
    mu.x[(9*i-2):(9*i)] <- ellipse$mu[1,i] + 2*pi
    
    mu.y[9*i - c(2,5,8)] <- ellipse$mu[2,i] - 2*pi
    mu.y[9*i - c(1,4,7)] <- ellipse$mu[2,i]
    mu.y[9*i - c(0,3,6)] <- ellipse$mu[2,i] + 2*pi
    
    d <- eigen(ellipse$cov.shape[[i]])
    a[(9*i-8):(9*i)] <- ellipse$r[i]*sqrt(d$values)[1]
    b[(9*i-8):(9*i)] <- ellipse$r[i]*sqrt(d$values)[2]
    angle[(9*i-8):(9*i)] <- atan(d$vectors[2,1]/d$vectors[1,1])
    
    
  }
  gg <- ggplot() +
    geom_ellipse(aes(x0 = mu.x, y0 = mu.y, a = a, b = b, angle = angle, fill=rep(color.array[1:k], rep(9,k))), alpha=0.5,  show.legend = F) +
    coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
    geom_point(aes(x=phi, y=psi), data_plot) +
    # ggtitle("Plot of ellipses") +
    # theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(~ phi)) + ylab(expression(~ psi))
  return(gg)
}

#visualization
visual.total <- function(data, conf, ncomp = 3, ellipse, alpha = 0.05){
  # Visualize conformal predictions set and ellipses
  # input : (n x 2 matrix) data, number of components K, (list) ellipse, level alpha
  # output : conformal predictions set and ellipses plot
  
  library(ggplot2)
  library(ggforce)
  data_plot <- as.data.frame(data)
  colnames(data_plot) <- c("phi", "psi")
  
  k <- length(ellipse$r)
  
  color.array <- c("yellow", "pink", "green", "cyan", "salmon", "brwon", "cornsilk", "olivegreen", "chocolate", "purple")
  mu.x <- matrix()
  mu.y <- c()
  a <- c()
  b <- c()
  angle <- c()
  
  for(i in 1:k){
    mu.x[(9*i-8):(9*i-6)] <- ellipse$mu[1,i] - 2*pi
    mu.x[(9*i-5):(9*i-3)] <- ellipse$mu[1,i]
    mu.x[(9*i-2):(9*i)] <- ellipse$mu[1,i] + 2*pi
    
    mu.y[9*i - c(2,5,8)] <- ellipse$mu[2,i] - 2*pi
    mu.y[9*i - c(1,4,7)] <- ellipse$mu[2,i]
    mu.y[9*i - c(0,3,6)] <- ellipse$mu[2,i] + 2*pi
    
    d <- eigen(ellipse$cov.shape[[i]])
    a[(9*i-8):(9*i)] <- ellipse$r[i]*sqrt(d$values)[1]
    b[(9*i-8):(9*i)] <- ellipse$r[i]*sqrt(d$values)[2]
    angle[(9*i-8):(9*i)] <- atan(d$vectors[2,1]/d$vectors[1,1])
    
    
  }
  gg <- ggplot() +
    annotate(geom = "raster", x = conf$phi, y = conf$psi,
             fill = scales::colour_ramp(c("green", "white"))(conf$pdf), interpolate = TRUE) +
    geom_ellipse(aes(x0 = mu.x, y0 = mu.y, a = a, b = b, angle = angle, fill=rep(color.array[1:k], rep(9,k))), alpha = 0.5,  show.legend = F) +
    coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
    geom_point(aes(x=phi, y=psi), data_plot) +
    # ggtitle(paste("vMM & ellipses (K=",  ncomp, ", alpha = ", alpha, ")", sep="")) +
    # theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(~ phi)) + ylab(expression(~ psi))
  return(gg)
}


######## vMM Method #####################################
vMM.method.parameter <- function(data, ncomp = 5){
  # Split data and fit von Mises mixture model
  # input : (n x 2 matrix) data, number of components K
  # output : (list) data, subdata for level, parameter, constant, number of components
  
  sp.data <- splitdata(data)
  parameter <- cal.parameter(sp.data$mix, ncomp)
  constant <- cal.constant(parameter)
  return(list(data = data, subdata = sp.data$level, parameter = parameter, constant = constant, ncomp = ncomp))
}

vMM.method.alpha <- function(first.cal, alpha = 0.05, grid.size = 100){
  # Obtain conformal prediction set for level alpha
  # input : list from vMM.method.parameter, level alpha, grid size
  # output : (list) data, subdata for level, parameter, constant, number of components
  
  level <- cal.level(first.cal$subdata, first.cal$parameter, first.cal$constant, alpha = alpha)
  level.max  <- cal.level(first.cal$subdata, first.cal$parameter, first.cal$constant, alpha = alpha, maximum = TRUE)
  conf <- cal.conformal.set(first.cal$parameter, first.cal$constant, level, grid.size)
  ellipse.parameter <- cal.ellipse(first.cal$parameter, first.cal$constant, level.max)
  return( list(data = first.cal$data,
               subdata = first.cal$subdata,
               parameter = first.cal$parameter,
               constant = first.cal$constant,
               ncomp = first.cal$ncomp,
               ellipse.parameter = ellipse.parameter, level = level, level.max = level.max, conformal.set = conf,
               vMM.plot = visual.vMM(first.cal$data, conf, first.cal$ncomp, alpha),
               ellipse.plot = visual.ellipse(first.cal$data, ellipse.parameter),
               total.plot = visual.total(first.cal$data, conf, first.cal$ncomp, ellipse.parameter, alpha)))
}

vMM.method <- function(data, alpha = 0.05, ncomp = 5,  grid.size = 100){
  # Obtain conformal prediction set for level alpha and number of components
  # input : (n x 2 matrix) data, level alpha, number of components, grid size
  # output : (list) results from vMM.method.parameter and vMM.method.alpha
  
  first.cal <- vMM.method.parameter(data, ncomp)
  return(vMM.method.alpha(first.cal, alpha))
}

######## Optimal K & Alpha #############################
vMM.optimization <- function(data, first.alpha = 0.1, grid.size = 100, alpha.size = 0.01){
  # Find optimal K and alpha and obtain conformal prediction set
  # input : (n x 2 matrix) data, initial alpha
  # output : (list) K plot, alpha plot, final conformal prediction set
  
  library(ggplot2)
  minimal.K_array <- c()
  minimal.alpha_array <- c()
  
  parameter.array <- list()
  constant.array <- list()
  
  sp.data <- splitdata(data)
  tryCatch(
    {
      for(K in 2:30){
        parameter.array[[K]] <- cal.parameter(sp.data$mix, K)
        constant.array[[K]] <- cal.constant(parameter.array[[K]])
        cat(paste("Fitting", K, "components is completed.\n"))
      }
    },
    error = function(cond) NA
  )
  
  
  K_array <- c()
  K.opt <- c()
  conf_array <- c()
  alpha.opt <- c()
  optimal.K <- c()
  minimum.K <- c()
  optimal.vMM.alpha <-c()
  minimum.vMM.alpha <-c()
  
  K.opt_array <- c()
  alpha.opt_array <- c()
  
  temp.alpha <- first.alpha
  gg.optimal.K <- c()
  gg.optimal.vMM.alpha <- c()
  alpha.opt <- c()
  optimal.vMM.method <- c()
  
  for(i in 1:10){
    K_array <- c()
    K.opt <- c()
    conf_array <- c()
    alpha.opt <- c()
    optimal.K <- c()
    minimum.K <- c()
    optimal.vMM.alpha <-c()
    minimum.vMM.alpha <-c()
    
    gg.optimal.K <- c()
    gg.optimal.vMM.alpha <- c()
    alpha.opt <- c()
    optimal.vMM.method <- c()
    
    for(K in 2:length(parameter.array)){
      level <- cal.level(sp.data$level, parameter.array[[K]], constant.array[[K]], alpha = temp.alpha)
      conf <- cal.conformal.set(parameter.array[[K]], constant.array[[K]], level, grid.size)
      K_array <- c(K_array, nrow(conf))
    }
    
    optimal.K <- cbind(2:(length(K_array)+1), K_array/grid.size^2)
    optimal.K <- as.data.frame(optimal.K)
    colnames(optimal.K) <- c("K", "set")
    minimum.K <- optimal.K[optimal.K$set==min(optimal.K$set),]
    if(length(minimum.K$K)>=2){
      K.opt <- minimum.K$K[1]
    } else{
      K.opt <- minimum.K$K
    }
    K.opt_array <- c(K.opt_array, K.opt)
    
    gg.optimal.K <- ggplot() +
      geom_line(aes(x=K, y=set), optimal.K, size=2) +
      geom_point(aes(x=K, y=set), minimum.K, color = 'blue', size=3) +
      # ggtitle("Plot of measure acoording to K") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab("K")
    
    tryCatch(
      {
        for(alpha in (1:floor(1/alpha.size/2))*alpha.size){
          conf_alpha <- vMM.method.alpha(list(data = data, subdata = sp.data$level,
                                              parameter = parameter.array[[K.opt]],
                                              constant = constant.array[[K.opt]],
                                              ncomp = K.opt),
                                         alpha = alpha, grid.size)
          conf_array <- c(conf_array, nrow(conf_alpha$conformal.set))
        }
      },
      error = function(cond) NA
    )
    
    optimal.vMM.alpha <- cbind((1:length(conf_array))*alpha.size, conf_array/grid.size^2)
    optimal.vMM.alpha <- as.data.frame(optimal.vMM.alpha)
    colnames(optimal.vMM.alpha) <- c("alpha", "set")
    minimum.vMM.alpha <- optimal.vMM.alpha[rowSums(optimal.vMM.alpha)==min(rowSums(optimal.vMM.alpha)),]
    if(length(minimum.vMM.alpha$alpha) >= 2){
      minimum.vMM.alpha <- minimum.vMM.alpha[1,]
    } else{
      minimum.vMM.alpha <- minimum.vMM.alpha
    }
    alpha.opt <- minimum.vMM.alpha$alpha
    alpha.opt_array <- c(alpha.opt_array, alpha.opt)
    
    gg.optimal.vMM.alpha <- ggplot() +
      geom_line(aes(x=alpha, y=set), optimal.vMM.alpha, size=2) +
      geom_abline(slope=-1, intercept=sum(minimum.vMM.alpha), size=2) +
      geom_point(aes(x=alpha, y=set), minimum.vMM.alpha, color = 'red', size=3) +
      coord_cartesian(xlim=c(0,0.5), ylim=c(0,0.5)) +
      # ggtitle("Plot of measure acoording to alpha") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab(expression(~ alpha))
    
    optimal.vMM.method <- vMM.method.alpha(list(data = data, subdata = sp.data$level,
                                                         parameter = parameter.array[[K.opt]],
                                                         constant = constant.array[[K.opt]],
                                                         ncomp = K.opt),
                                                    alpha = alpha.opt)
    
    if(temp.alpha == alpha.opt){
      break
    } else{
      temp.alpha <- alpha.opt
    }
  }
  
  optimal.K.alpha <- as.data.frame(matrix(c(K.opt, alpha.opt),1,2))
  colnames(optimal.K.alpha) <- c("K", "alpha")

  return(list(optimal.K.alpha = optimal.K.alpha,
              final = optimal.vMM.method,
              # parameter.array = parameter.array,
              # constant.array = constant.array,
              # K.opt_array = K.opt_array,
              # alpha.opt_array = alpha.opt_array,
              K.plot = gg.optimal.K,
              alpha.plot = gg.optimal.vMM.alpha))
}

####### Clustering ##################################
Clustering.vMM <- function(optimization){
  # Clustering by von Mises mixture method
  # input : (n x 2 matrix)  data, (n/2 x 2 matrix) subdata for level, (6 x K matrix) parameter
  # output : clustering plots
   
  data <- optimization$final$data
  subdata <- optimization$final$subdata
  parameter <- optimization$final$parameter
  alpha <- optimization$optimal.K.alpha$alpha
  
  Conf.Scores <- conf.score.bvn.approx(subdata, parameter)
  Conf.Pred.Set <- conf.pred.set(subdata, Conf.Scores$conf.score,
                                     parameter, alpha = alpha, grid.size = 200)
  Index.Clusters <- conf.pred.clusters(parameter = parameter, conf.score.level = Conf.Pred.Set$Level)
  Cluster.Assignment <- conf.pred.assign.clusters(x = data, parameter = parameter,
                                                           conf.score.level = Conf.Pred.Set$Level,
                                                           Index.Clusters = Index.Clusters,
                                                           option.clusters = 1)
  Cluster.Assignment.op2 <- conf.pred.assign.clusters(x = data, parameter = parameter,
                                                               conf.score.level = Conf.Pred.Set$Level,
                                                               Index.Clusters = Index.Clusters,
                                                               option.clusters = 2)
  plot1 <- Draw.ellipses.bvn.approx(data = data, parameter.ellipse = parameter,
                                    conf.score.level = Conf.Pred.Set$Level,
                                    cut.point = Conf.Pred.Set$Point,
                                    assigned.clusters = Cluster.Assignment$Cluster.Assignment[,3])
  plot2 <- Draw.ellipses.bvn.approx(data = data, parameter.ellipse = parameter,
                                    conf.score.level = Conf.Pred.Set$Level,
                                    cut.point = Conf.Pred.Set$Point,
                                    assigned.clusters = Cluster.Assignment.op2$Cluster.Assignment[,3])
  return(list(all.plot = plot1, zero.plot = plot2))
}


