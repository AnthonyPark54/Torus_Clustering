
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
# grid.size = 100
# first.alpha = 0.1
###########################################################

### basic functions ###

cal.ker <- function(x, y, kap1 = 25, kap2 = 25){
  # Calculate kernel for each point
  # input : (2 vector) a point, (n x 2 matrix) data, kappa 1, kappa 2
  # output : a kernel
  
  ker1 <- exp(kap1*cos(x[1]-y[,1]))
  ker2 <- exp(kap2*cos(x[2]-y[,2]))
  return(sum(ker1*ker2)/nrow(y)/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0)))
}

cal.ker.vector <- function(x, y, kap1 = 25, kap2 = 25){
  # Calculate kernel for set
  # input : (n x 2 vector) point set, (n x 2 matrix) data, kappa 1, kappa 2
  # output : kernel of set
  
  n <- nrow(x)
  ker.vector <- rep(0, n)
  for(i in 1:n){
    ker.vector[i] <- cal.ker(x[i,], y, kap1, kap2)
  }
  return(ker.vector)
}

# L.minus, plus
cal.level.mp <- function(x, alpha = 0.05, kap1 = 25, kap2 = 25){
  # Calculate the level of L-, L+
  # input : (n x 2 matrix) data, level alpha, kappa 1, kappa 2
  # output : level of L-, L+
  
  n <- nrow(x)
  ialpha <- floor((n+1)*alpha)
  den_est <- cal.ker.vector(x, x, kap1, kap2)
  ord_den_est <- den_est[order(den_est)]
  
  level_minus <- ord_den_est[ialpha]
  level_plus <- ord_den_est[ialpha] - (exp(kap1)*exp(kap2)-exp(-kap1)*exp(-kap2))/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0))/n
  return(c(level_minus, level_plus))
}

conformal.L.mp <- function(x, alpha = 0.05, kap1 = 25, kap2 =25, grid.size = 100){
  # Obtain L-, L+
  # input : (n x 2 matrix) data, level alpha, kappa 1, kappa 2, grid size
  # output : L-, L+. between set
  
  level <- cal.level.mp(x, alpha, kap1, kap2)
  names <- c("phi","psi", "pdf")
  Axis <- seq(0, 2*pi, length = grid.size)
  lattice <- cbind(rep(Axis, grid.size),rep(Axis, each=grid.size))
  
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
  
  return(list(data = x, minus = L_minus, plus = L_plus, btw = btw.L,
              alpha = alpha, kap1 = kap1, kap2 = kap2, grid.size = grid.size))
}

# conformal set
cal.sigma.level <- function(x, lattice, kap1 = 25, kap2 = 25){
  # Calculate the level of conformal prediction set with KDE
  # input : (n x 2 matrix) data, level alpha, kappa 1, kappa 2
  # output : ;evel of conformal prediction set with KDE
  
  n <- nrow(x)
  ker.pdf <- cal.ker.vector(lattice, x, kap1, kap2)
  return((ker.pdf*n+exp(kap1+kap2)/(2*pi*besselI(kap1,0))/(2*pi*besselI(kap2,0)))/(n+1))
}

conformal.prediction <- function(L.mp){
  # Obtain conformal prediction set with KDE
  # input : (list) L.mp from conformal L.mp
  # output : conformal prediction set with KDE
  
  names <- c("phi", "psi", "pdf")
  n <- nrow(L.mp$data)
  m <- nrow(L.mp$btw)
  
  lattice <- L.mp$btw[,1:2]
  sigma.level <- cal.sigma.level(L.mp$data, lattice, L.mp$kap1, L.mp$kap2)
  
  pi.y <- rep(0, m)
  for(i in 1:m){
    new.x <- rbind(L.mp$data,lattice[i,])
    pi.y[i] <- (sum(cal.ker.vector(L.mp$data, new.x, L.mp$kap1, L.mp$kap2) <= sigma.level[i])+1)/(n+1)
  }
  
  C.alpha <- cbind(lattice[which(pi.y > L.mp$alpha),], L.mp$btw[,3][which(pi.y > L.mp$alpha)])
  C.alpha <- rbind(as.matrix(L.mp$minus), C.alpha)
  C.alpha <- as.data.frame(C.alpha)
  colnames(C.alpha) <- names
  
  return(list (L.mp = L.mp, C.alpha = C.alpha))
}

# visualization
visual.L <- function(L.mp, only.plus = FALSE, only.minus = FALSE){
  # Visualize L-, L+
  # input : (list) L.mp from conformal L.mp, only plus, only minus
  # output : L-, L+ plots
  
  library(ggplot2)
  data_plot <- as.data.frame(L.mp$data)
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

visual.conformal <- function(conf, with.plus = FALSE, with.minus = FALSE, only.conf = FALSE){
  # Visualize L-, L+, conformal prediction set
  # input : (list) L.mp from conformal.prediction, with plus, with minus, only conf
  # output : L-, L+, conformal prediction set plots
  
  library(ggplot2)
  data_plot <- as.data.frame(conf$L.mp$data)
  colnames(data_plot) <- c("phi", "psi")
  
  if(with.plus){
    gg <- ggplot() +
      annotate(geom = "raster", x = conf$L.mp$plus$phi, y = conf$L.mp$plus$psi,
               fill = scales::colour_ramp(c("green", "white"))(conf$L.mp$plus$pdf), interpolate = TRUE) +
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
      annotate(geom = "raster", x = conf$C.alpha$phi, y = conf$C.alpha$psi,
               fill = scales::colour_ramp(c("green", "white"))(conf$C.alpha$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) + 
      geom_raster(aes(x=phi, y=psi, fill=pdf), conf$L.mp$minus, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set with L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(only.conf){
    gg <- ggplot(conf$C.alpha, aes(x=phi, y=psi)) + geom_raster(aes(fill=pdf), interpolate = TRUE) +
      scale_fill_gradientn(colours=c("green","white"), guide = FALSE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  
  if(!with.plus & !with.minus & !only.conf){
    gg <- ggplot() +
      annotate(geom = "raster", x = conf$L.mp$plus$phi, y = conf$L.mp$plus$psi,
               fill = scales::colour_ramp(c("green", "white"))(conf$L.mp$plus$pdf), interpolate = TRUE) +
      coord_cartesian(xlim=c(0,2*pi), ylim=c(0,2*pi)) +
      annotate(geom = "raster", x = conf$C.alpha$phi, y = conf$C.alpha$psi,
               fill = scales::colour_ramp(c("purple", "white"))(conf$C.alpha$pdf), interpolate = TRUE) +
      geom_raster(aes(x=phi, y=psi, fill=pdf), conf$L.mp$minus, interpolate = TRUE) +
      scale_fill_gradientn(colours=c("pink","white"), guide = FALSE) +
      geom_point(aes(x=phi, y=psi), data_plot) +
      # ggtitle("Plot of conformal prediction set with L plus & L minus") +
      # theme(plot.title = element_text(hjust = 0.5)) +
      xlab(expression(~ phi)) + ylab(expression(~ psi))
  }
  return(gg)
}


######## Kernel Method #####################################
kernel.method <- function(data, alpha = 0.05, kap1 = 25, kap2 = 25,  grid.size = 100, conformal.set = FALSE){
  # Obtain conformal prediction set with kernel
  # input : (n x 2 matrix) torus data, level alpha, kappa 1, kappa 2, grid.size, with conformal set
  # output : (list) conformal prediction set and plots
  
  L <- conformal.L.mp(data, alpha, kap1, kap2, grid.size)
  if(conformal.set){
    conf <- conformal.prediction(L)
    return(list(L.set = L, L.plot = visual.L(L), conformal.set = conf$C.alpha, conformal.plot = visual.conformal(conf)))
  } else if(!conformal.set){
    return(list(L.set = L, L.plot = visual.L(L)))
  }
}

######## Optimal Kappa & Alpha #############################
kernel.optimization <- function(data, first.alpha = 0.1, grid.size = 100, alpha.size = 0.01, kappa.size = 1){
  # Find optimal kappa, alpha
  # input : (n x 2 matrix) torus data, intial alpha, grid size
  # output : (list) optimal conformal set
  
  Lplus_array <- c()
  Lminus_array <- c()
  conf_array <- c()
  kappa_array <- (1:floor(100/kappa.size))*kappa.size
  
  kappa.opt <- c()
  alpha.opt <- c()
  
  temp.alpha <- first.alpha
  gg.optimal.kappa <- c()
  gg.optimal.kernel.alpha <- c()
  optimal.kernel.method <- c()
  
  kappa.opt_array <- c()
  alpha.opt_array <- c()
  
  for(i in 1:10){
    Lplus_array <- c()
    Lminus_array <- c()
    conf_array <- c()
    kappa.opt <- c()
    alpha.opt <- c()
    gg.optimal.kappa <- c()
    gg.optimal.kernel.alpha <- c()
    optimal.kernel.method <- c()
    
    for(kappa in kappa_array){
      conf_kappa <- kernel.method(data, conformal.set = TRUE, alpha = first.alpha,
                                  kap1 = kappa, kap2 = kappa, grid.size = grid.size)
      Lplus_array <- c(Lplus_array, nrow(conf_kappa$L.set$plus))
      Lminus_array <- c(Lminus_array, nrow(conf_kappa$L.set$minus))
      conf_array <- c(conf_array, nrow(conf_kappa$conformal.set))
      cat(paste("Conformal prediction with kappa ", kappa, " is completed.\n"))
    }
    optimal.kappa <- cbind(kappa_array, Lplus_array/grid.size^2, conf_array/grid.size^2, Lminus_array/grid.size^2)
    optimal.kappa <- as.data.frame(optimal.kappa)
    colnames(optimal.kappa) <- c("kappa", "Lplus", "conf", "Lminus")
    minimum.kappa <- optimal.kappa[optimal.kappa$conf==min(optimal.kappa$conf),]
    if(length(minimum.kappa$kappa)>=2){
      kappa.opt <- minimum.kappa$kappa[1]
    } else{
      kappa.opt <- minimum.kappa$kappa
    }
    kappa.opt_array <- c(kappa.opt_array, kappa.opt)
    
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
      {for(alpha in (1:floor(0.3/alpha.size))*alpha.size){
        conf_alpha <- kernel.method(data, conformal.set = TRUE, alpha = alpha, kap1 = kappa.opt, kap2 = kappa.opt, grid.size = grid.size)
        minus_array <- c(minus_array, nrow(conf_alpha$L.set$minus))
        plus_array <- c(plus_array, nrow(conf_alpha$L.set$plus))
        conf_array <- c(conf_array, nrow(conf_alpha$conformal.set))
        cat(paste("Conformal prediction with ", alpha, " is completed.\n"))
      }},
      error = function(cond) NA
    )
    
    optimal.kernel.alpha <- cbind((1:length(minus_array))*alpha.size, plus_array/grid.size^2, conf_array/grid.size^2, minus_array/grid.size^2)
    optimal.kernel.alpha <- as.data.frame(optimal.kernel.alpha)
    colnames(optimal.kernel.alpha) <- c("alpha", "Lplus", "conf", "Lminus")
    minimum.kernel.alpha <- optimal.kernel.alpha[rowSums(optimal.kernel.alpha[,c(1,3)])==min(rowSums(optimal.kernel.alpha[,c(1,3)])),]
    alpha.opt <- minimum.kernel.alpha$alpha
    alpha.opt_array <- c(alpha.opt_array, alpha.opt)
    
    gg.optimal.kernel.alpha <- ggplot() +
      geom_line(aes(x=alpha, y=conf,  colour = "Conformal"), optimal.kernel.alpha, size = 2) +
      geom_line(aes(x=alpha, y=Lplus,  colour = 'Lplus'), optimal.kernel.alpha, size = 1) +
      geom_line(aes(x=alpha, y=Lminus,  colour = 'Lminus'), optimal.kernel.alpha, size = 1) +
      coord_cartesian(xlim=c(0,0.3), ylim=c(0,0.3)) +
      scale_color_discrete(name = "Set") +
      theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
      geom_abline(slope=-1, intercept=sum(minimum.kernel.alpha[c(1,3)]), size =1) +
      geom_point(aes(x=alpha, y=conf), minimum.kernel.alpha, color = 'black', size=3) +
      ylab(expression(tilde(mu)(widehat(C)^(alpha)))) + xlab(expression(~ alpha))
    
    optimal.kernel.method <- kernel.method(data, conformal.set = TRUE,  alpha = alpha.opt,
                                                    kap1 = kappa.opt, kap2 = kappa.opt, grid.size = grid.size)
    
    if(temp.alpha == alpha.opt){
      break
    } else{
      temp.alpha <- alpha.opt
    }
  }
  
  optimal.kappa.alpha <- as.data.frame(matrix(c(kappa.opt, alpha.opt),1,2))
  colnames(optimal.kappa.alpha) <- c("kappa", "alpha")
  
  return(list(optimal.kappa = optimal.kappa,
              minimum.kappa = minimum.kappa,
              kappa.plot = gg.optimal.kappa,
              kappa.opt_array = kappa.opt_array,
              alpha.opt_array = alpha.opt_array,
              optimal.kernel.alpha = optimal.kernel.alpha,
              minimum.kernel.alpha = minimum.kernel.alpha,
              alpha.plot = gg.optimal.kernel.alpha,
              final = optimal.kernel.method))
  
}






