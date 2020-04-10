
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
# div.num = 100
# first.alpha = 0.1
###########################################################

### basic functions ###

splitdata <- function(x){
  n <- nrow(x)
  index <- sample(n, floor(n/2))
  x1 <- x[index,]
  x2 <- x[-index,]
  
  return(list(mix = x1, level = x2))
}
cal.parameter <- function(x1, ncomp = 3){
  library(BAMBI)
  summ <- summary(fit_vmsinmix(x1, method = "rwmh", ncomp = ncomp, unimodal.component = TRUE, n.iter = 1000))
  return(summ$`estimate `$mode)
}
cal.constant <- function(parameter){
  k <- ncol(parameter)
  d <- parameter
  
  div.num <- 200
  x_axis <- seq(0, 2*pi, length=div.num)
  lattice <- cbind(rep(x_axis,div.num),rep(x_axis,each=div.num))
  
  const <- rep(0, k)
  
  y <- matrix(0, nrow(lattice), k)
  for(i in 1:k){
    for(j in 1:nrow(lattice)){
      y[j,i] <- exp(d[2,i]*cos(lattice[j,1]-d[5,i])
                    + d[3,i]*cos(lattice[j,2]-d[6,i])
                    + d[4,i]*sin(lattice[j,1]-d[5,i])*sin(lattice[j,2]-d[6,i]))
    }
    const[i] <- div.num^2/sum(y[,i])/4/pi/pi
  }
  return(const)
}
cal.pdf.each <- function(lattice, parameter, constant){
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
cal.level <- function(x2, parameter, constant, alpha = 0.05, maximum = FALSE){
  pdf.each <- cal.pdf.each(x2, parameter, constant)
  if(maximum){
    pdf <- apply(pdf.each, 1, max)
  } else{
    pdf <- rowSums(pdf.each)
  }
  
  n <- nrow(x2)
  ialpha <- floor((n+1)*alpha)
  return(pdf[order(pdf)][ialpha])
}

# VMM
cal.conformal.set <- function(parameter, constant, pdf.level, div.num = 100){
  k <- ncol(parameter)
  d <- parameter
  names <- c("phi","psi", "pdf")
  
  x_axis <- seq(0, 2*pi, length=div.num)
  lattice <- cbind(rep(x_axis,div.num),rep(x_axis,each=div.num))
  
  y <- cal.pdf.each(lattice, parameter, constant)
  pdf <- rowSums(y)
  conf <- cbind(lattice[which(pdf >= pdf.level),], pdf[which(pdf >= pdf.level)])
  conf <- as.data.frame(conf)
  colnames(conf) <- names
  return(conf)
}
visual.vMM <- function(data, conf, ncomp = 3, alpha = 0.1){
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
  sp.data <- splitdata(data)
  parameter <- cal.parameter(sp.data$mix, ncomp)
  constant <- cal.constant(parameter)
  return(list(data = data, subdata = sp.data$level, parameter = parameter, constant = constant, ncomp = ncomp))
}
vMM.method.alpha <- function(first.cal, alpha = 0.05, div.num = 100){
  level <- cal.level(first.cal$subdata, first.cal$parameter, first.cal$constant, alpha = alpha)
  level.max  <- cal.level(first.cal$subdata, first.cal$parameter, first.cal$constant, alpha = alpha, maximum = TRUE)
  conf <- cal.conformal.set(first.cal$parameter, first.cal$constant, level, div.num)
  ellipse.parameter <- cal.ellipse(first.cal$parameter, first.cal$constant, level.max)
  return( list(ellipse.parameter = ellipse.parameter, level = level, level.max = level.max, conformal.set = conf,
               vMM.plot = visual.vMM(first.cal$data, conf, first.cal$ncomp, alpha),
               ellipse.plot = visual.ellipse(first.cal$data, ellipse.parameter),
               total.plot = visual.total(first.cal$data, conf, first.cal$ncomp, ellipse.parameter, alpha)))
}

vMM.method <- function(data, alpha = 0.05, ncomp = 5,  div.num = 100){
  first.cal <- vMM.method.parameter(data, ncomp)
  second.cal <- vMM.method.alpha(first.cal, alpha)
  return(list(first.cal, second.cal))
}

######## Optimal K & Alpha #############################
vMM.optimization <- function(data, first.alpha = 0.1){
  library(ggplot2)
  K_array <- c()
  sp.data <- splitdata(data)
  tryCatch(
    {
      for(K in 2:30){
        parameter <- cal.parameter(sp.data$mix, K)
        constant <- cal.constant(parameter)
        level <- cal.level(sp.data$level, parameter, constant, alpha = first.alpha)
        conf <- cal.conformal.set(parameter, constant, level, div.num = 100)
        K_array <- c(K_array, nrow(conf))
      }
    },
    error = function(cond) NA
  )
  
  optimal.K <- cbind(2:(length(K_array)+1), K_array/100^2)
  optimal.K <- as.data.frame(optimal.K)
  colnames(optimal.K) <- c("K", "set")
  minimum.K <- optimal.K[optimal.K$set==min(optimal.K$set),]
  if(length(minimum.K$K)>=2){
    K.opt <- minimum.K$K[1]
  } else{
    K.opt <- minimum.K$K
  }
  
  gg.optimal.K <- ggplot() +
    geom_line(aes(x=K, y=set), optimal.K, size=2) +
    geom_point(aes(x=K, y=set), minimum.K, color = 'blue', size=3) +
    # ggtitle("Plot of measure acoording to K") +
    # theme(plot.title = element_text(hjust = 0.5)) +
    ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab("K")

  conf_array <- c()
  oror.parm <- vMM.method.parameter(data, ncomp = K.opt)
  tryCatch(
    {
      for(alpha in (1:25)/50){
        oror <- vMM.method.alpha(oror.parm, alpha = alpha)
        conf_array <- c(conf_array, nrow(oror$conformal.set))
      }
    },
    error = function(cond) NA
  )
  
  optimal.vMM.alpha <- cbind((1:length(conf_array))/50, conf_array/100^2)
  optimal.vMM.alpha <- as.data.frame(optimal.vMM.alpha)
  colnames(optimal.vMM.alpha) <- c("alpha", "set")
  minimum.vMM.alpha <- optimal.vMM.alpha[rowSums(optimal.vMM.alpha)==min(rowSums(optimal.vMM.alpha)),]
  if(length(minimum.vMM.alpha$alpha)>=2){
    minimum.vMM.alpha <- minimum.vMM.alpha[1,]
  } else{
    minimum.vMM.alpha <- minimum.vMM.alpha
  }
  vMM.alpha.opt <- minimum.vMM.alpha$alpha
  
  
  gg.optimal.vMM.alpha <- ggplot() +
    geom_line(aes(x=alpha, y=set), optimal.vMM.alpha, size=2) +
    geom_abline(slope=-1, intercept=sum(minimum.vMM.alpha), size=2) +
    geom_point(aes(x=alpha, y=set), minimum.vMM.alpha, color = 'red', size=3) +
    coord_cartesian(xlim=c(0,0.5), ylim=c(0,0.5)) +
    # ggtitle("Plot of measure acoording to alpha") +
    # theme(plot.title = element_text(hjust = 0.5)) +
    ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab(expression(~ alpha))
  
  original.optimal.vMM.method <- vMM.method.alpha(oror.parm, alpha = vMM.alpha.opt)
  original.optimal.vMM.method$vMM.plot
  
  return(list(optimal.K = optimal.K,
              minimum.K = minimum.K,
              K.plot = gg.optimal.K,
              optimal.vMM.alpha = optimal.vMM.alpha,
              minimum.vMM.alpha = minimum.vMM.alpha,
              vMM.alpha.opt = vMM.alpha.opt,
              alpha.plot = gg.optimal.vMM.alpha,
              final = original.optimal.vMM.method,
              optimal.plot = original.optimal.vMM.method$vMM.plot))
}
