
# This function conducts conformal prediction sets on the torus data with kernel density method.
# Finally, conduct L-, L+, C_alpha on the torus data.
# For clustering, this function find optimal prediction set.

###########################################################
## Torus Data : n x 2 dataframe ###########################
#### each row : [0,2*pi) x [0,2*pi) #######################
###########################################################
## Parameter setting ######################################
# alpha = 0.05
# kap1 = 25, kap2 = 25
# div.num = 100
# first.alpha = 0.1
###########################################################

### basic functions ###

cal.ker <- function(x, y, kap1 = 25, kap2 = 25){
  ker1 <- exp(kap1*cos(x[1]-y[,1]))
  ker2 <- exp(kap2*cos(x[2]-y[,2]))
  return(sum(ker1*ker2)/nrow(y)/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0)))
}
cal.ker.vector <- function(x, y, kap1 = 25, kap2 = 25){
  n <- nrow(x)
  ker.vector <- rep(0, n)
  for(i in 1:n){
    ker.vector[i] <- cal.ker(x[i,], y, kap1, kap2)
  }
  return(ker.vector)
}

# L.minus, plus
cal.level.mp <- function(x, alpha = 0.05, kap1 = 25, kap2 = 25){
  n <- nrow(x)
  ialpha <- floor((n+1)*alpha)
  den_est <- cal.ker.vector(x, x, kap1, kap2)
  ord_den_est <- den_est[order(den_est)]
  
  level_minus <- ord_den_est[ialpha]
  level_plus <- ord_den_est[ialpha] - (exp(kap1)*exp(kap2)-exp(-kap1)*exp(-kap2))/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0))/n
  return(c(level_minus, level_plus))
}
conformal.L.mp <- function(x, alpha = 0.05, kap1 = 25, kap2 =25, div.num = 100){
  level <- cal.level.mp(x, alpha, kap1, kap2)
  names <- c("phi","psi", "pdf")
  x_axis <- seq(0, 2*pi, length=div.num)
  lattice <- cbind(rep(x_axis,div.num),rep(x_axis,each=div.num))
  
  ker.pdf <- cal.ker.vector(lattice, x, kap1, kap2)
  
  L_minus <- cbind(lattice[which(ker.pdf >= level[1]),], ker.pdf[which(ker.pdf >= level[1])])
  L_minus <- as.data.frame(L_minus)
  colnames(L_minus) <- names

  L_plus <- cbind(lattice[which(ker.pdf >= level[2]),], ker.pdf[which(ker.pdf >= level[2])])
  L_plus <- as.data.frame(L_plus)
  colnames(L_plus) <- names
  
  btw.L <- cbind(lattice[which((ker.pdf >= level[2]) & (ker.pdf < level[1])),],
                 ker.pdf[which((ker.pdf >= level[2]) & (ker.pdf < level[1]))])
  colnames(btw.L) <- names
  
  return(list(minus = L_minus, plus = L_plus, btw = btw.L))
}

# conformal set
cal.sigma.level <- function(x, lattice, kap1 = 25, kap2 = 25, div.num = 100){
  x_axis <- seq(0, 2*pi, length=div.num)
  n <- nrow(x)
  ker.pdf <- cal.ker.vector(lattice, x, kap1, kap2)
  return((ker.pdf*n+exp(kap1+kap2)/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0)))/(n+1))
}
conformal.prediction <- function(x, L.mp, alpha = 0.05, kap1 = 25, kap2 = 25, div.num = 100){
  
  names <- c("phi", "psi", "pdf")
  n <- nrow(x)
  m <- nrow(L.mp$btw)
  
  lattice <- L.mp$btw[,1:2]
  sigma.level <- cal.sigma.level(x, lattice, kap1, kap2, div.num)
  
  pi.y <- rep(0, m)
  for(i in 1:m){
    new.x <- rbind(x,lattice[i,])
    pi.y[i] <- (sum(cal.ker.vector(x ,new.x, kap1, kap2) <= sigma.level[i])+1)/(n+1)
  }
  
  C.alpha <- cbind(lattice[which(pi.y > alpha),], L.mp$btw[,3][which(pi.y > alpha)])
  C.alpha <- rbind(as.matrix(L.mp$minus), C.alpha)
  C.alpha <- as.data.frame(C.alpha)
  colnames(C.alpha) <- names
  
  return(C.alpha)
}

# visualization
visual.L <- function(data, L.mp, only.plus = FALSE, only.minus = FALSE){
  library(ggplot2)
  data_plot <- as.data.frame(data)
  colnames(data_plot) <- c("phi", "psi")
  
  if(only.plus){
    gg <- ggplot(L.mp$plus, aes(x=phi, y=psi)) +  geom_raster(aes(fill=pdf), interpolate = TRUE) +
      scale_fill_gradientn(colours=c("green","white"), guide = FALSE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of L plus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(only.minus){
    gg <- ggplot(L.mp$minus, aes(x=phi, y=psi)) + geom_raster(aes(fill=pdf), interpolate = TRUE) + 
      scale_fill_gradientn(colours=c("green","white"), guide = FALSE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(!only.plus & !only.minus){
    gg <- ggplot() +
      annotate(geom = "raster", x = L.mp$plus$phi, y = L.mp$plus$psi,
               fill = scales::colour_ramp(c("green", "white"))(L.mp$plus$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_raster(aes(x=phi, y=psi, fill=pdf), L.mp$minus, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of L plus & L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  return(gg)
}
visual.conformal <- function(data, conf, L.mp, with.plus = FALSE, with.minus = FALSE, only.conf = FALSE){
  library(ggplot2)
  data_plot <- as.data.frame(data)
  colnames(data_plot) <- c("phi", "psi")
  
  if(with.plus){
    gg <- ggplot() +
      annotate(geom = "raster", x = L.mp$plus$phi, y = L.mp$plus$psi,
               fill = scales::colour_ramp(c("green", "white"))(L.mp$plus$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_raster(aes(x=phi, y=psi, fill=pdf), conf, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set with L plus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(with.minus){
    gg <- ggplot() +
      annotate(geom = "raster", x = conf$phi, y = conf$psi,
               fill = scales::colour_ramp(c("green", "white"))(conf$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) + 
      geom_raster(aes(x=phi, y=psi, fill=pdf), L.mp$minus, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set with L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(only.conf){
    gg <- ggplot(conf, aes(x=phi, y=psi)) + geom_raster(aes(fill=pdf), interpolate = TRUE) +
      scale_fill_gradientn(colours=c("green","white"), guide = FALSE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(!with.plus & !with.minus & !only.conf){
    gg <- ggplot() +
      annotate(geom = "raster", x = L.mp$plus$phi, y = L.mp$plus$psi,
               fill = scales::colour_ramp(c("green", "white"))(L.mp$plus$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      annotate(geom = "raster", x = conf$phi, y = conf$psi,
               fill = scales::colour_ramp(c("purple", "white"))(conf$pdf), interpolate = TRUE) +
      geom_raster(aes(x=phi, y=psi, fill=pdf), L.mp$minus, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set with L plus & L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  return(gg)
}


######## Kernel Method #####################################
kernel.method <- function(data, conformal.set = FALSE, alpha = 0.05, kap1 = 25, kap2 = 25,  div.num = 100){
  L <- conformal.L.mp(data, alpha, kap1, kap2, div.num)
  if(conformal.set){
    conf <- conformal.prediction(data, L, alpha, kap1, kap2, div.num)
    return(list(L.set = L, L.plot = visual.L(data, L), conformal.set = conf, conformal.plot = visual.conformal(data, conf, L)))
  } else if(!conformal.set){
    return(list(L.set = L, L.plot = visual.L(data, L)))
  }
}

######## Optimal Kappa & Alpha #############################
kernel.optimization <- function(data, first.alpha = 0.1){
  Lplus_array <- c()
  Lminus_array <- c()
  conf_array <- c()
  kappa_array <- 1:100
  for(kappa in kappa_array){
    oror <- kernel.method(data, conformal.set = TRUE, alpha = first.alpha, kap1 = kappa, kap2 = kappa, div.num = 100)
    Lplus_array <- c(Lplus_array, nrow(oror$L.set$plus))
    Lminus_array <- c(Lminus_array, nrow(oror$L.set$minus))
    conf_array <- c(conf_array, nrow(oror$conformal.set))
  }
  optimal.kappa <- cbind(kappa_array, Lplus_array/100^2, conf_array/100^2, Lminus_array/100^2)
  optimal.kappa <- as.data.frame(optimal.kappa)
  colnames(optimal.kappa) <- c("kappa", "Lplus", "conf", "Lminus")
  minimum.kappa <- optimal.kappa[optimal.kappa$conf==min(optimal.kappa$conf),]
  if(length(minimum.kappa$kappa)>=2){
    kappa.opt <- minimum.kappa$kappa[1]
  } else{
    kappa.opt <- minimum.kappa$kappa
  }
  
  
  gg.optimal.kappa <- ggplot() +
    geom_line(aes(x=kappa, y=conf,  color = "Conformal"), optimal.kappa, size=2) +
    geom_line(aes(x=kappa, y=Lplus, color="Lplus"), optimal.kappa, size=1) +
    geom_line(aes(x=kappa, y=Lminus,  color = "Lminus"), optimal.kappa, size=1) +
    scale_color_discrete(name = "Set") +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
    geom_point(aes(x=kappa, y=conf), minimum.kappa, color = 'black', size=3) +
    ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab(expression(~ kappa))
  
  minus_array <- c()
  plus_array <- c()
  conf_array <- c()
  tryCatch(
    {for(alpha in (1:50)/100){
        oror <- kernel.method(data, conformal.set = TRUE, alpha = alpha, kap1 = kappa.opt, kap2 = kappa.opt, div.num = 100)
        minus_array <- c(minus_array, nrow(oror$L.set$minus))
        plus_array <- c(plus_array, nrow(oror$L.set$plus))
        conf_array <- c(conf_array, nrow(oror$conformal.set))
    }},
    error = function(cond) NA
  )
  
  optimal.kernel.alpha <- cbind((1:length(minus_array))/100, plus_array/100^2, conf_array/100^2, minus_array/100^2)
  optimal.kernel.alpha <- as.data.frame(optimal.kernel.alpha)
  colnames(optimal.kernel.alpha) <- c("alpha", "Lplus", "conf", "Lminus")
  minimum.kernel.alpha <- optimal.kernel.alpha[rowSums(optimal.kernel.alpha[,c(1,3)])==min(rowSums(optimal.kernel.alpha[,c(1,3)])),]
  kernel.alpha.opt <- minimum.kernel.alpha$alpha
  
  
  gg.optimal.kernel.alpha <- ggplot() +
    geom_line(aes(x=alpha, y=conf,  colour = "Conformal"), optimal.kernel.alpha, size = 2) +
    geom_line(aes(x=alpha, y=Lplus,  colour = 'Lplus'), optimal.kernel.alpha, size = 1) +
    geom_line(aes(x=alpha, y=Lminus,  colour = 'Lminus'), optimal.kernel.alpha, size = 1) +
    coord_cartesian(xlim=c(0,0.5), ylim=c(0,0.5)) +
    scale_color_discrete(name = "Set") +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
    geom_abline(slope=-1, intercept=sum(minimum.kernel.alpha[c(1,3)]), size =1) +
    geom_point(aes(x=alpha, y=conf), minimum.kernel.alpha, color = 'black', size=3) +
    ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab(expression(~ alpha))
  
  original.optimal.kernel.method <- kernel.method(data, conformal.set = TRUE,  alpha = kernel.alpha.opt,
                                                  kap1 = kappa.opt, kap2 = kappa.opt, div.num = 100)
  return(list(optimal.kappa = optimal.kappa,
              minimum.kappa = minimum.kappa,
              kappa.plot = gg.optimal.kappa,
              optimal.kernel.alpha = optimal.kernel.alpha,
              minimum.kernel.alpha = minimum.kernel.alpha,
              alpha.plot = gg.optimal.kernel.alpha,
              final = original.optimal.kernel.method,
              optimal.plot = original.optimal.kernel.method$conformal.plot))
  
}






