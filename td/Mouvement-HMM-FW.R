# HMM pour mouvement animal

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(mclust)

################################################################################
# Données
data <- read.table('dryad_zebra.csv', header=TRUE, sep=',')
position <- as.matrix(data[, c('x', 'y')])
Y <- as.matrix(data[, c('speed', 'angle')]); n <- nrow(Y)
Y <- scale(Y)

# Parms
mu = list(-0.5*c(1, 1), c(0, 0), 0.5*c(1, 1))
K <- length(mu)
sigma <- list(0.5*diag(2), 0.5*diag(2), 0.5*diag(2))
nu <- rep(1/K, K)
pi <- matrix(0.1, K, K); diag(pi) <- 1 - rowSums(pi)
parms = list(nu=nu, pi=pi, mu=mu, sigma=sigma)
parms

################################################################################
# Fonctions HMM
Forward <- function(Y, parms){
  K <- length(parms$mu); n <- nrow(Y)
  # Densités d'émission
  phi <- matrix(0, n, K)
  for(k in 1:K){phi[, k] <- dmvnorm(Y, mean=parms$mu[[k]], sigma=parms$sigma[[k]], log=TRUE)}
  phi[which(phi < -100)] <- -100; phi <- exp(phi)
  # Récurrence avant
  Fw <- matrix(0, n, K); logL <- 0
  Fw[1,] <- phi[1,] * parms$nu
  logL <- logL + log(sum(Fw[1,]))
  Fw[1,] <- Fw[1,] / sum(Fw[1,])
  for (t in (2:n)){
    Fw[t,] <- Fw[t-1,] %*% parms$pi
    Fw[t,] <- Fw[t,] * phi[t,]
    logL <- logL + log(sum(Fw[t,]))
    Fw[t,] <- Fw[t,] / sum(Fw[t,])
  }
  return(list(phi=phi, Fw=Fw, logL=logL))
}

Backward <- function(Y, parms, forward){
  K <- length(parms$mu); n <- nrow(Y)
  # Récurence arrière
  tau <- G <- matrix(0, n, K)
  tau[n,] <- forward$Fw[n, ]
  # eStep$eta(k, l) <- sum_t E(Z(t, k) Z(t+1, l) | X)
  eta <- matrix(0, K, K)
  for (t in ((n-1):1)){
    G[t+1, ] <- forward$Fw[t,] %*% parms$pi
    B <- as.vector(tau[t+1, ] / G[t+1, ])
    tau[t,] <- forward$Fw[t,] * (parms$pi %*% B)
    tau[t,] <- tau[t,] / sum(tau[t,])
    eta <- eta + parms$pi * outer(forward$Fw[t,], B)
  }
  return(list(tau=tau, eta=eta))
}

forward <- Forward(Y, parms)
backward <- Backward(Y, parms, forward)
