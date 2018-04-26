###########################
#### Full conditionals ####
###########################

#l es el limite ese de la serie geometrica
#### Acceptance probability

#j para el peso. CHecar funcion y ver si guardar los vjs. Los pesos son secuenciales, ie, hazlo funcion y velo corriendo
acc_prob <- function(param, x, tau, mu, wj, j, new_vj, old_vj, new_tau, new_mu, l){
  switch(param,
         vj = {
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu)^2)
           expo_arg <- exp(-(tau/2)*expo_arg)
           expo <- expo_arg * wj
           denom <- sum(log(sum(sapply(1-rowSums(expo), function(x) x^(1:l))) + 1))
           
           old_vj[j] <- new_vj
           new_w <- c(old_vj[1],  old_vj[-1]*cumprod(1-old_vj[-length(old_vj)]))
           expo <- expo_arg * new_w
           num <- sum(log(sum(sapply(1 - rowSums(expo), function(x) x^(1:l))) + 1))
           prob <- min(1, exp(num-denom))
           runif(1) < prob
         },
         tau = {
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu)^2)
           num <- exp(-(new_tau/2)*expo_arg)
           num <- num * wj
           denom <- exp(-(tau/2)*expo_arg)
           denom <- denom * wj
           denom <- sum(log(sum(sapply(1 - rowSums(denom), function(x) x^(1:l))) + 1))
           num <- sum(log(sum(sapply(1 - rowSums(num), function(x) x^(1:l)) ) + 1))
           prob <- min(1, exp(num-denom))
           runif(1) < prob
         },
         mu = {
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu[-j])^2)
           expo_arg <- exp(-(tau/2)*expo_arg)
           expo <- expo_arg * wj[-j]
           expo_num <- expo + (exp(-(tau/2)*(x-new_mu)^2)*wj[j])
           expo_denom <- expo + (exp(-(tau/2)*(x-mu[j])^2)*wj[j]) 
           num <- sum(log(sum(sapply(1 - rowSums(expo_num), function(x) x^(1:l)) ) + 1))
           denom <- sum(log(sum(sapply(1 - rowSums(expo_denom), function(x) x^(1:l)) ) + 1))
          
           prob <- min(1, exp(num-denom))
           runif(1) < prob
         })
}

acc_prob <- cmpfun(acc_prob)
#### w_j

v_j <- function(j, b, d){
  nj <- length(which(d == j))
  nj_p <- length(which(d > j))
  rbeta(1, 1 + nj, b + nj_p)
}

v_j <- cmpfun(v_j)

wj <- function(x, di, b, m, tau, mu, old_vjs, old_wj, l){
  vjs_prop <- sapply(1:m, v_j, b = b, d = di)
  for(i in 1:length(vjs_prop)){
    si <- acc_prob('vj', x, tau, mu, j=i, wj=old_wj,
                   new_vj = vjs_prop[i], old_vj = old_vjs, l = l)
    if(si) 
      old_vjs[i] <- vjs_prop[i]
  }
  wjs <- vector('double', length(old_vjs))
  wjs[1] <- old_vjs[1]
  wjs[-1] <- old_vjs[-1]*cumprod(1 - old_vjs[-length(old_vjs)])
  # if(any(is.na(wjs)))
  #   message('Hay NAs en los pesos')
  list(wj=wjs, vj=old_vjs)
}

wj <- cmpfun(wj)

#### \rho_j 

rho_j <- function(x, old_rho, R, d_i, mu, tau, m){
 cl <- makeCluster(4)
 registerDoParallel(cl)
 probs <- foreach(i = 1:m, .combine = 'cbind', .export = 'dens') %dopar% {
   vapply(R, dens, 0, j = i, tau=tau, x=x, mu=mu, d_i=d_i)
 }
 stopCluster(cl) 
  probs <- apply(probs, 2, function(x){
    sorted <- sort(x, index.return = T)
    indices <- sorted$ix
    sorted <- sorted$x
    sorted <- exp(sorted - sorted[length(sorted)])
    norms <- sorted/sum(sorted)
    norms[indices]
  })
  
  apply(probs, 2, function(x) sample(R, 1, prob = x))
}

rho_j <- cmpfun(rho_j)

dens <- function(r, j, tau, x, mu, d_i){
  sumable <- which(d_i == j)
  suma <- sum(x[sumable + 1] - mu[j] - r*(x[sumable] - mu[j])^2/(1-r^2))
  log(1-r^2)*(-length(sumable)/2) + -(tau/2)*suma
}

dens <- cmpfun(dens)
#### d_i


d_i <- function(x, old_d, wj, mu, tau, rho){
  sigma <- sqrt(1/tau); n <- length(x) - 1
  dens <- function(i, j){ 
    log(wj[j]) - (x[i+1] - mu[j])^2/(2*sigma^2) - (x[i] -  (mu[j]+rho[j]*(x[i+1]-mu[j]))^2/(2*sigma^2))
  }

  probs <- matrix(nrow = length(wj), ncol = n)
  for(i in 1:n){
    probs[,i] <- vapply(1:length(wj), dens, 0, i = i)
  }
  probs <- apply(probs, 2, function(x){
    sorted <- sort(x, index.return = T)
    indices <- sorted$ix
    sorted <- sorted$x
    sorted <- exp(sorted - sorted[length(sorted)])
    norms <- sorted/sum(sorted)
    norms[indices]
  })
  apply(probs, 2, function(x) sample(1:length(wj), 1, prob = x))
}

d_i <- cmpfun(d_i)

#### z_li

z_li <- function(x, old_z, wj, mu, sigma, ki, di, error){
 
  dens <- function(i, j){ 
    log(wj[j]) + log(1-exp(-(x[i]-mu[j])^2/(2*sigma)))
  }

  probs <- matrix(nrow = length(wj), ncol = length(x)-1)
  for(i in 1:ncol(probs)){
    probs[,i] <- vapply(1:length(wj), dens, 0, i = i) 
  }
  probs <- apply(probs, 2, function(x){
    sorted <- sort(x, index.return = T)
    indices <- sorted$ix
    sorted <- sorted$x
    sorted <- sorted - sorted[length(sorted)]
    sorted <- ifelse(sorted >= log(error) - log(length(sorted)), exp(sorted), 0)
    norms <- sorted/sum(sorted)
    norms[indices]
  })
  lapply(1:ncol(probs), function(i) sample(1:length(wj), ki[i], T, probs[,i]))
}


#### mu_j

mu_j <- function(x, old_mu, rho, d_i, t, m.hyper, tau, wj, l){
  medias_prop <- sapply(sort(unique(d_i)), mu.prop, d_i=d_i, tau=tau, rho=rho,t=t,m=m.hyper, x=x)
  for(i in 1:length(medias_prop)){
    si <- acc_prob('mu', x, tau, mu = old_mu, j=i, wj=wj,
                   new_mu = medias_prop[i],l=l)
    if(si) 
      old_mu[i] <- medias_prop[i]
  }
  old_mu
}
mu_j <- cmpfun(mu_j)

mu.prop <- function(j, d_i, rho, tau,m,t, x){

  indices <- which(d_i == j)
  frac <- tau/(1-rho[j])
  nj <- length(indices)
  suma <- sum((x[indices + 1] + x[indices])/(1+rho[j]))
  varianza <- t + (frac*nj) + (nj*tau)
  media <- ((m*t) + (frac*suma)) / varianza
  rnorm(1, media, sqrt(varianza))
}
mu.prop <- cmpfun(mu.prop)

#### tau

tau <- function(x, mu, d_i, old_tau, rho, a, b, wj, l){
  n <- length(x)-1
  ambos <- x[-n] - mu[d_i]
  a.hat <- a + n/2
  suma <- (x[-1] - mu[d_i] - (rho[d_i]*ambos))^2/(1-rho[d_i]^2)
  b.hat <- b + ambos^2/2 + suma/2
  tau_prop <- rgamma(1, a.hat, rate = b.hat)
  si <- acc_prob('tau', x, old_tau, mu = mu, wj=wj,
                 new_tau = tau_prop, l=l)
  if(1/tau_prop == 0){
    stop('Underflow en la varianza')
  }
  
  if(si)
    return(tau_prop)
  else
    return(old_tau)
}

tau <- cmpfun(tau)
### k_i

k_i <- function(old_k, x, wj, p, z_il, mu, sigma){
  z.new <- sample(1:length(wj), length(old_k), T, prob = wj)
  comp <- vector('double', length(x) - 1)
  for(i in 1:length(comp)){
    comp[i] <- x[i+1] - mu[z.new[i]]
  }
  prob <- p/(1-p) * (1 - exp(-1/(2*sigma) * comp^2) )
  prob <- sapply(prob, function(x) min(1, x))
  ratios <- runif(length(prob))
  
  comp.des <- vector('double', length(x) - 1)
  for(i in 1:(length(x)-1)){
    comp.des[i] <- x[i+1] - mu[  z_il[[i]][old_k[i]]  ]
  }
  prob.des <- (1-p)/p * 1/(1 - exp(-1/(2*sigma) * comp.des^2) )
  prob.des <- sapply(prob.des, function(x) min(1, x))
  
  no.acepto <- sapply(1:length(prob), function(i) ratios[i] > prob[i])
  new.k <- numeric(length(old_k))
  for(i in 1:length(no.acepto)){
    if(no.acepto[i] && old_k[i] > 1){
        r <- runif(1)
        if(r < prob.des[i])
          new.k[i] <- old_k[i] - 1
        else
          new.k[i] <- old_k[i]
    }
    else{
      new.k[i] <- old_k[i] + 1
      z_il[[i]] <- append(z_il[[i]], z.new[i])
    }
  } 
  eval(parse(text=paste("z <- list(", paste(z_il, collapse=','), ')')),
             envir = parent.frame())# :O
  return(new.k)
}

# rtruncgamma <- function(shape, rate, prev.tau, trunc){
#   trunc <- mpfr(trunc, prec = 100); rate <- mpfr(rate, prec=100)
#   y <- runif(1)*exp(-rate*mpfr(prev.tau, 100))
#   if(y == 0)
#     y <- mpfr(1e-320, 100)
#   x <- mpfr(runif(1), prec = 100)
#   (x * (trunc^shape + (-log(y))^shape )/ (shape*(shape-1))  )^(1/shape)
# }
##########


# simulacion <- function(n){
#   dens <- function(x) 0.2*dnorm(x, -5, 1)+0.4*dnorm(x, 0, 3)+0.4*dnorm(x, 7) 
#   dens <- Vectorize(dens, 'x')
#   curve(dens, xlim=c(-15,15))
#   
#   componentes <- sample(1:3,prob=c(0.2, 0.4, 0.4),size=n,replace=TRUE)
#   mus <- c(-5, 0, 7)
#   sds <- sqrt(c(1, 3, 1))
#   
#   samples <- rnorm(n, mean=mus[componentes], sd=sds[componentes])
#   return(samples)
# }

#ademas de tau o sigma
mcmc <- function(datos, nsim,  R, b, c, a, t, m, m.hyper, init, l){
  library(foreach); library(doParallel)
  w <- vector('list', nsim); w[[1]] <- init[['w']]
  old_vjs <- vector('list', nsim); old_vjs[[1]] <- init[['vjs']]
  d <- init[['d']]
  #z <- init[['z']]
  rho <- vector('list', nsim); rho[[1]] <- init[['rho']]
  mu <- vector('list', nsim); mu[[1]] <- init[['mu']]
  tau <- vector('list', nsim); tau[[1]] <- init[['sigma']]
#  k <- vector('list', nsim); k[[1]] <- init[['k']]
 pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    #list(R, mu, dj, z_il, sigma)
  for(i in 2:nsim){
    #print(i)
   if(i %% 100 == 0)
     setTxtProgressBar(pb, i)
    
    pesos <- wj(x = datos, d, b, m, tau[[i-1]], mu[[i-1]], 
                old_vjs[[i-1]],old_wj=w[[i-1]], l=l)
    w[[i]] <- pesos[['wj']]
    old_vjs[[i]] <- pesos[['vj']]
    
    d <- d_i(x = datos, old_d = d, 
            wj = w[[i]], mu = mu[[i-1]],
            tau = tau[[i-1]], rho = rho[[i-1]])
    rho[[i]] <- rho_j(x = datos, old_rho = rho[[i-1]], R,
                d_i = d, mu = mu[[i-1]],
                tau = tau[[i-1]], m=m)
    mu[[i]] <- mu_j(x = datos, old_mu = mu[[i-1]],
              rho = rho[[i]], d_i = d, m = m, t =t, wj = w[[i]],
              tau = tau[[i-1]], l=l)
    tau[[i]] <- tau(x = datos, mu = mu[[i]], d_i = d,
          old_tau = tau[[i-1]], a=a, b=b, wj=w[[i]],rho = rho[[i]], l=l)
  }
 close(pb)
  return(list(w, rho, mu, tau))
}
mcmc <- cmpfun(mcmc)


post.inf <- function(mcmc, burn_in){
  w.hat <- colMeans(do.call(cbind, mcmc[['w']])[,-(1:burnin)])
  m.hat <- colMeans(do.call(cbind, mcmc[['mu']])[,-(1:burnin)])
  sigma.hat <- colMeans(do.call(cbind, mcmc[['sigma']])[,-(1:burnin)])
  rho.hat <- colMeans(do.call(cbind, mcmc[['rho']])[,-(1:burnin)])
  return(list(w.hat, m.hat, sigma.hat, rho.hat))
}

VJS <- rbeta(5,1,1)
W <- 1:5; W[1] <- VJS[1]; W[2:5] <- VJS[-1]*cumprod(1-VJS[-5])
init <- list(d = sample(1:5, length(datos$promedios)-1, T),
             w = W,
             vjs = VJS,
             rho = rep(0.5, 5),
             sigma = 2,
             mu = rep(1, 5))  
#nota: checar x_0
itera <- mcmc(datos$promedios, 20000, 
              R = seq(0.001, 0.999, by = 0.001),
              a = 5,l=10 , m = 5, b = 1, init = init, m.hyper= 300,
              t = 1) 



# rmodelo <- function(x, rho, sigma, mus, wj){
#   pesos <- wj*dnorm(mpfr(sim_x[1], 100), mus, sqrt(sigma))
#   pesos <- asNumeric(pesos/sum(pesos))
#   cual <- sample(1:length(wj), 1, prob = pesos)
#   rnorm(1, mus[cual] + rho*(x - mus[cual]), sqrt((1-rho[cual]^2)*sigma))
# }
# 
# sim_x <- x_n
# for(i in 2:length(sim_x)){
#   if(i %% 100 == 0){
#     print(i)
#   }
#   sim_x[i] <- rmodelo(sim_x[i-1], sigma = varianza,
#                       rho = rhos, mus = medias, wj = weights)
# }



