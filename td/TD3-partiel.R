# TD3 : Réseau trophique

rm(list=ls())

# Data
Y <- read.table('chilean_TI.csv', header=TRUE, sep=';')
rownames(Y) <- Y[, 2]; 
Y <- as.matrix(Y[, -(1:2)])
diag(Y) <- 0
n <- nrow(Y)
image(1:n, 1:n, Y)
K <- 5
data <- list(Y=Y)

# Init
library(mclust)
init <- Mclust(Y, G=K)
tau <- init$z
veStep <- list(tau=tau)

VEstep <- function(data, mStep, veStep){
  n <- nrow(data$Y); K <- length(mStep$pi)
  logPhi <- array(dim=c(n, n, K, K))
  for(k in 1:K){for( l in 1:K){
    logPhi[, , k, l] <- dbinom(data$Y, 1, prob=mStep$gamma[k, l], log=TRUE)
    diag(logPhi[, , k, l]) <- 0
  }}
  tau <- matrix(0, n, K)
  for(i in 1:n){for (k in 1:K){
    tau[i, k] <- mStep$pi[k]*exp(sum(logPhi[i, , k, ]*veStep$tau))
  }}
  return(list(tau=tau, logPhi=logPhi))
}
Mstep <- function(data, veStep){
  pi <- colMeans(veStep$tau)
  J_I <- matrix(1, nrow(data$Y), ncol(data$Y)) - diag(nrow(data$Y))
  gamma <- t(veStep$tau)%*%data$Y%*%veStep$tau / t(veStep$tau)%*%J_I%*%veStep$tau
  return(list(pi=pi, gamma=gamma))
}

mStep <- Mstep(data=data, veStep=veStep)
veStep <- VEstep(data=data, mStep=mStep, veStep=veStep)
