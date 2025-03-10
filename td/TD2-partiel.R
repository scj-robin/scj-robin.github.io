# TD2 : Mouvement animal

rm(list=ls())
library(mvtnorm)

# Data 
data <- read.table('dryad_zebra.csv', heade=TRUE, sep=',')
head(data)
Y <- as.matrix(data[, c('speed', 'angle')])
X <- as.matrix(data[, c('x', 'y')])
n <- nrow(Y)
K <- 3
data <- list(X=X, Y=Y)
plot(X, type='b', pch=20)

# Init
library(mclust)
init <- Mclust(Y, G=K)
tau <- init$z
eta <- array(0, dim=c(n, K, K))
for(t in 2:n){eta[t, , ] <- tau[t-1, ]%o%tau[t, ]}
eStep <- list(tau=tau, eta=eta)

Forward <- function(data, mStep){
  n <- nrow(data$Y); K <- length(mStep$nu)
  phi <- t(sapply(1:n, function(t){sapply(1:K, function(k){
    mvtnorm::dmvnorm(data$Y[t, ], mean=mStep$mu[k, ], sigma=mStep$sigma[k, , ])
  })}))
  F <- matrix(NA, n, K)
  F[1, ] <- mStep$nu * phi[1, ]
  logLik <- log(sum(F[1, ]))
  F[1, ] <- F[1, ] / sum(F[1, ])
  for(t in 2:n){
    F[t, ] <- phi[t, ] * (F[t-1, ]%*%mStep$pi)
    logLik <- logLik + log(sum(F[t, ]))
    F[t, ] <- F[t, ] / sum(F[t, ])
  }
  return(list(F=F, phi=phi, logLik=logLik))
}
Backward <- function(forward, mStep){
  n <- nrow(forward$F); K <- length(mStep$nu)
  tau <- matrix(NA, n, K); eta <- array(0, dim=c(n, K, K))
  tau[n, ] <- forward$F[n, ]
  for(t in (n-1):1){
    G <- as.vector(forward$F[t, ]%*%mStep$pi)
    B <- tau[t+1, ]/G
    eta[t+1, , ] <- (forward$F[t, ]%o%B) * mStep$pi
    tau[t, ] <- rowSums(eta[t+1, , ])
    tau[t, ] <- tau[t, ] / sum(tau[t, ])
  }
  return(list(tau=tau, eta=eta))
}
Estep <- function(data, mStep){
  forward <- Forward(data=data, mStep=mStep)
  backward <- Backward(forward=forward, mStep=mStep)
  return(list(tau=backward$tau, eta=backward$eta, phi=forward$phi, logLik=forward$logLik))
}
Mstep <- function(data, eStep){
  nu <- eStep$tau[1, ]
  pi <- apply(eStep$eta, c(2, 3), sum)
  pi <- pi / rowSums(pi)
  N <- colSums(eStep$tau); K <- length(N)
  mu <- t(eStep$tau)%*%data$Y / N%o%rep(1, 2)
  sigma <- array(NA, dim=c(K, ncol(Y), ncol(Y)))
  for(k in 1:K){
    sigma[k, , ] <- t(data$Y)%*%diag(eStep$tau[, k])%*%data$Y / N[k]
    sigma[k, , ] <- sigma[k, , ] - mu[k, ]%o%mu[k, ]
  }
  return(list(nu=nu, pi=pi, mu=mu, sigma=sigma))
}

mStep <- Mstep(data=data, eStep=eStep)
forward <- Forward(data=data, mStep=mStep)
