###########################
#### Full conditionals ####
###########################

#l es el limite ese de la serie geometrica
#### Acceptance probability

#j para el peso. CHecar funcion y ver si guardar los vjs. Los pesos son secuenciales, ie, hazlo funcion y velo corriendo
acc_prob <- function(param, x, tau, mu, wj, j, new_vj, old_vj, new_tau, new_mu, l, b, d_j, a.hat, b.hat, b.beta, media_full, tau_full){
  switch(param,
         vj = {
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu)^2)
           expo_arg <- exp(-(tau/2)*expo_arg)
           expo <- expo_arg * wj
           denom <- sum(log(sapply(1-colSums(expo), function(x) sum(x^(1:l))) + 1))
           
           #Lo que faltaba
           nj <- length(which(d_j == j))
           nj_p <- length(which(d_j > j))
           falta <- 0
           if(dbeta(new_vj, 1 + nj, b.beta + nj_p) == 0){
             if(dbeta(old_vj[j], 1 + nj, b.beta + nj_p) == 0)
               falta <- 0
             else
               falta <- -log(dbeta(old_vj[j], 1 + nj, b.beta + nj_p))
           }
           falta <- falta - dbeta(new_vj, 9*old_vj[j], 9*(1-old_vj[j]), log=T) + dbeta(old_vj[j], 9*new_vj, 9*(1-new_vj),log=T)
           
           old_vj[j] <- new_vj
           new_w <- c(old_vj[1],  old_vj[-1]*cumprod(1-old_vj[-length(old_vj)]))
           expo <- expo_arg * new_w
           num <- sum(log(sapply(1 - colSums(expo), function(x) sum(x^(1:l))) + 1))
           stopifnot(length(num) == 1)
           
           prob <- min(1, exp(num-denom + falta))
           runif(1) < prob
         },
         tau = {
           falta <- 0
           if(dgamma(new_tau, a.hat, rate = b.hat) == 0){
             if(dgamma(tau, a.hat, rate = b.hat) == 0)
               falta <- 0
             else
               falta <- -dgamma(tau, a.hat, rate = b.hat, log = T)
           } else{
             if(dgamma(tau, a.hat, rate = b.hat) == 0)
               falta <- dgamma(new_tau, a.hat, rate = b.hat, log = T)
             else
               falta <- dgamma(new_tau, a.hat, rate = b.hat, log = T) - dgamma(tau, a.hat, rate = b.hat, log = T)
           }
             
           falta <- falta + dgamma(tau, 90*new_tau, rate = 90, log = T) - dgamma(new_tau, 90*tau, rate = 90, log= T)
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu)^2)
           num <- exp(-(new_tau/2)*expo_arg)
           num <- num * wj
           denom <- exp(-(tau/2)*expo_arg)
           denom <- denom * wj
           denom <- sum(log(sapply(1 - colSums(denom), function(x) sum(x^(1:l))) + 1))
           num <- sum(log(sapply(1 - colSums(num), function(x) sum(x^(1:l)) ) + 1))
           prob <- min(1, exp(num-denom+falta))
           runif(1) < prob
         },
         mu = {
           m <- length(wj)
           x <- x[-length(x)]
           expo_arg <- sapply(x, function(y) (y-mu[-j])^2)
           expo_arg <- exp(-(tau/2)*expo_arg)
           
           falta <- -(tau_full/2)*((new_mu-media_full)^2 - (mu[j]-media_full)^2)
           
           expo <- expo_arg * wj[-j]
           expo_num <- expo + (exp(-(tau/2)*(x-new_mu)^2)*wj[j])
           expo_denom <- expo + (exp(-(tau/2)*(x-mu[j])^2)*wj[j]) 
           num <- sum(log(sapply(1 - colSums(expo_num), function(x) sum(x^(1:l)) ) + 1))
           denom <- sum(log(sapply(1 - colSums(expo_denom), function(x) sum(x^(1:l)) ) + 1))
           if(is.nan(num) || is.nan(denom)){
             num <- prod(sapply(1 - colSums(expo_num), function(x) sum(x^(1:l)) ) + 1)
             denom <- prod(sapply(1 - colSums(expo_denom), function(x) sum(x^(1:l)) ) + 1)
             if(num == 0)
               prob <- 0
             else if(denom == 0)
               prob <- 1
             else
               prob <- min(1, exp(falta)*num/denom)
           } else
                prob <- min(1, exp(num-denom+falta))
           runif(1) < prob
         })
}

acc_prob <- cmpfun(acc_prob)
#### w_j

v_j <- function(j, b, d){
  nj <- length(which(d == j))
  nj_p <- length(which(d > j))
  dbeta(x, 1 + nj, b + nj_p)
}

v_j <- cmpfun(v_j)

wj <- function(x, di, b.beta, m, tau, mu, old_vjs, old_wj, l, empty){
  vjs_prop <- sapply(old_vjs, function(x) rbeta(1,9*x, 9*(1-x)))
#  vjs_prop <- sapply(1:m, v_j, b = b, d = di)
  for(i in 1:length(vjs_prop)){
    if(empty[i]){
      old_vjs[i] <- rbeta(1, 1, b.beta)
      next
    }
    si <- acc_prob('vj', x, tau, mu, j=i, wj=old_wj,
                   new_vj = vjs_prop[i], old_vj = old_vjs, l = l, b.beta = b.beta, d_j = di)
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

rho_j <- function(x, old_rho, R, d_i, mu, tau, m, empty){
  cl <- makeCluster(2)
  registerDoParallel(cl)
  probs <- foreach(i = 1:m, .combine = 'cbind', .export = 'dens') %dopar% {
    if(empty[i])
      rep(1/length(R),length(R))
    else
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
  rho_post <- apply(probs, 2, function(x) sample(R, 1, prob = x))
  rho_post[empty] <- sample(R, sum(empty))
  rho_post
}

rho_j <- cmpfun(rho_j)

dens <- function(r, j, tau, x, mu, d_i){
  sumable <- which(d_i == j)
  suma <- sum((x[sumable + 1] - mu[j] - r*(x[sumable] - mu[j]))^2/(1-r^2))
  (log(1-r^2)*(-length(sumable)/2)) -((tau/2)*suma)
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
   new_d <- apply(probs, 2, function(x) sample(1:length(wj), 1, prob = x))
   empty <- sapply(1:length(wj), function(i) !any(new_d == i))
   list(new_d = new_d, empty = empty)
   #stopifnot(length(unique(new_d)) == 3)
   #return(new_d)
}

d_i <- cmpfun(d_i)

#### z_li



#### mu_j

mu_j <- function(x, old_mu, rho, d_i, t, m.hyper, tau, wj, l, empty){
  par_full <- sapply(1:length(wj), mu.prop, d_i=d_i, tau=tau, rho=rho,t=t,m=m.hyper, x=x, empty=empty)
  g <- function(new_mu, a.hat, b.hat, x, wj, l, mu, j, old_mu){
    old_mu[j] <- new_mu

    falta <- dnorm(new_mu, a.hat, 1/sqrt(b.hat), log = T)
    
    x <- x[-length(x)]
    expo_arg <- sapply(x, function(y) (y-old_mu)^2)
    num <- exp(-(mu/2)*expo_arg)
    num <- num * wj
    num_b <- sum(log( sapply(1 - colSums(num), function(x) sum(x^(1:l)) ) + 1))
    if(is.nan(num_b))
      num_b <- log(prod(sapply(1 - colSums(num), function(x) sum(x^(1:l)) ) + 1))
    num_b + falta
  }
  g <- cmpfun(g)
  
  for(i in 1:length(old_mu)){
    if(empty[i]){
      old_mu[i] <- rnorm(1, m.hyper, 1/sqrt(t))
      next
    }
    
    old_mu[i] <- uni.slice(old_mu[i], g, a.hat = par_full[1,i], b.hat = par_full[2,i], mu = tau, x = x, wj = wj, l=l, j = i, old_mu = old_mu)
  }
    # si <- acc_prob('mu', x, tau, mu = old_mu, j=i, wj=wj, media_full = par_full[1,i],
    #                tau_full = par_full[2,i], new_mu = medias_prop[i],l=l)
    # if(si){
    #   if(i == 1)
    #     eval(parse(text = 'contador.mu <- contador.mu + 1'), envir=.GlobalEnv)
    #   old_mu[i] <- medias_prop[i] 
    # }
     
  
  old_mu
}
mu_j <- cmpfun(mu_j)

mu.prop <- function(j, d_i, rho, tau,m,t, x, empty){
  if(empty[j])
    return(c(NA, NA))
  indices <- which(d_i == j)
  frac <- tau/(1-rho[j])
  nj <- length(indices)
  suma <- sum((x[indices + 1] + x[indices])/(1+rho[j]))
  varianza <- t + (frac*nj) + (nj*tau)
  media <- ((m*t) + (frac*suma)) / varianza
  c(media, varianza)
}
mu.prop <- cmpfun(mu.prop)

#### tau

# tau <- function(x, mu, d_i, old_tau, rho, a, b, wj, l){
#   n <- length(x)-1
#   ambos <- x[-n] - mu[d_i]
#   a.hat <- a + n/2
#   suma <- (x[-1] - mu[d_i] - (rho[d_i]*ambos))^2/(1-rho[d_i]^2)
#   b.hat <- b + sum(ambos)^2/2 + sum(suma)/2
#   tau_prop <- rgamma(1, 100*old_tau, rate = 100)
#   #tau_prop <- rgamma(1, a.hat, rate = b.hat)
#   si <- acc_prob('tau', x, old_tau, mu = mu, wj=wj,
#                  new_tau = tau_prop, l=l, a.hat=a.hat, b.hat=b.hat)
#   if(1/tau_prop == 0){
#     stop('Underflow en la varianza')
#   }
#   
#   if(si){
#     eval(parse(text = 'contador.tau <- contador.tau + 1'), envir=.GlobalEnv)
#     return(tau_prop)
#   }
#   else
#     return(old_tau)
# }

tau <- function(x, mu, d_i, old_tau, rho, a, b, wj, l){
  n <- length(x)-1
  ambos <- x[-n] - mu[d_i]
  a.hat <- a + n/2
  suma <- (x[-1] - mu[d_i] - (rho[d_i]*ambos))^2/(1-rho[d_i]^2)
  b.hat <- b + sum(ambos)^2/2 + sum(suma)/2
  g <- function(new_tau, a.hat, b.hat, x, wj, l, mu, j, old_mu){
    falta <- dgamma(new_tau, a.hat, rate = b.hat, log = T)
    x <- x[-length(x)]
    expo_arg <- sapply(x, function(y) (y-mu)^2)
    num <- exp(-(new_tau/2)*expo_arg)
    num <- num * wj
    num_b <- sum(log(sapply(1 - colSums(num), function(x) sum(x^(1:l)) ) + 1))
    if(is.nan(num_b))
      num_b <- log(prod(sapply(1 - colSums(num), function(x) sum(x^(1:l)) ) + 1))
    num_b + falta
  }
  g <- cmpfun(g)
  
  uni.slice(old_tau, g, lower = 0, a.hat = a.hat, b.hat = b.hat, mu = mu, x = x, wj = wj, l=l)
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
mcmc <- function(datos, nsim,  R, b, b.beta, c, a, t, m, m.hyper,init, l){
  library(foreach); library(doParallel)
  w <- vector('list', nsim); w[[1]] <- init[['w']]
  old_vjs <- vector('list', nsim); old_vjs[[1]] <- init[['vjs']]
  
  
  d <- vector('list', nsim); d[[1]] <- init[['d']]
  #z <- init[['z']]
  rho <- vector('list', nsim); rho[[1]] <- init[['rho']]
  mu <- vector('list', nsim); mu[[1]] <- init[['mu']]
  tau <- vector('list', nsim); tau[[1]] <- init[['sigma']]
  empty <- matrix(ncol = nsim, nrow = length(w[[1]]))
  empty[,1] <- rep(F,nrow(empty))
#  k <- vector('list', nsim); k[[1]] <- init[['k']]
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    #list(R, mu, dj, z_il, sigma)
  for(i in 2:nsim){
    #print(i)
   if(i %% 100 == 0)
    setTxtProgressBar(pb, i)
    
    d_ <- d_i(x = datos, old_d = d[[i-1]], 
             wj = w[[i-1]], mu = mu[[i-1]],
             tau = tau[[i-1]], rho = rho[[i-1]])
    empty[,i] <- d_[['empty']]
    d[[i]] <- d_[['new_d']]
    
    pesos <- wj(x = datos, d[[i]], b.beta, m, tau[[i-1]], mu[[i-1]], 
                old_vjs[[i-1]],old_wj=w[[i-1]], l=l, empty = empty[,i])
    w[[i]] <- pesos[['wj']]
    old_vjs[[i]] <- pesos[['vj']]
    
    rho[[i]] <- rho_j(x = datos, old_rho = rho[[i-1]], R,
                d_i = d[[i]], mu = mu[[i-1]],
                tau = tau[[i-1]], m=m, empty = empty[,i])
    mu[[i]] <- mu_j(x = datos, old_mu = mu[[i-1]],
              rho = rho[[i]], d_i = d[[i]], m.hyper = m.hyper, t =t, wj = w[[i]],
              tau = tau[[i-1]], l=l, empty = empty[,i])
    tau[[i]] <- tau(x = datos, mu = mu[[i]], d_i = d[[i]],
          old_tau = tau[[i-1]], a=a, b=b, wj=w[[i]],rho = rho[[i]], l=l)
  }
  close(pb)
  return(list(w, rho, mu, tau, empty, d))
}
mcmc <- cmpfun(mcmc)


# post.inf <- function(mcmc, burnin){
#   w.hat <- rowMeans(do.call(cbind, mcmc[[1]])[,-(1:burnin)])
#   m.hat <- rowMeans(do.call(cbind, mcmc[[3]])[,-(1:burnin)])
#   sigma.hat <- mean(unlist(mcmc[[4]])[-(1:burnin)])
#   rho.hat <- rowMeans(do.call(cbind, mcmc[[2]])[,-(1:burnin)])
#   return(list(w.hat, m.hat, sigma.hat, rho.hat))
# }

post.inf <- function(mcmc, burnin){
  w.hat <- do.call(cbind, mcmc[[1]])
  m.hat <- do.call(cbind, mcmc[[3]])
  sigma.hat <- unlist(mcmc[[4]])
  rho.hat <- do.call(cbind, mcmc[[2]])
  return(list(w.hat, m.hat, sigma.hat, rho.hat))
}

#grupo <- kmeans(datos$promedios,5)
VJS <- rbeta(15,1,1)
W <- 1:15; W[1] <- VJS[1]; W[2:15] <- VJS[-1]*cumprod(1-VJS[-15])
init <- list(d = kmeans(o3[-1],15)$cluster,
             w = W,
             vjs = VJS,
             rho = rep(0.5, 15),
             sigma = 0.001,
             mu = kmeans(o3[-1],15)$centers)  

#pi <- read.csv('~/Desktop/probas.csv2')[,2]
#nota: checar x_0
#162
itera <- mcmc(o3, 100000, 
              R = seq(0.001, 0.999, by = 0.001),
              a = 0.9,l=10 , m = 15, b =  1, b.beta = 0.7, init = init, m.hyper= 60,
              t = 0.001111111) 


library(Rmpfr)
rmodelo <- function(x, rho, tau, mus, wj){
  sigma <- 1/tau
  pesos <- wj*dnorm(mpfr(x,100), mus, sqrt(sigma))
  pesos <- pesos/sum(pesos)
  res <- sum(pesos * mus)
  asNumeric(res)
}


#### Label switching

estimaciones <- post.inf(itera, 69999)

mcmc.pars <- array(data = NA, dim = c(1000, 15, 3) )
mcmc.pars[,,1] <- t(estimaciones[[2]])[90001:100000,][(1:10000)%%10 == 0,] #mu
mcmc.pars[,,2] <- t(estimaciones[[1]])[90001:100000,][(1:10000)%%10 == 0,]#wj
mcmc.pars[,,3] <- t(estimaciones[[4]])[90001:100000,][(1:10000)%%10 == 0,] #rho

  
 
complete.normal.loglikelihood <- function(x, z, pars){
  mu <- pars[,1]
  wj <- pars[,2]
  d_i <- z
  rho <- pars[,3]
  
  x_i <- x
  x_menos <- c(57, x[-length(x)])
  
  expo_arg <- sapply(x_menos, function(y) (y-mu)^2)
  denom <- exp(-( 0.652384/2)*expo_arg)
  denom <- denom * wj
  denom <- sum(log(sapply(1 - colSums(denom), function(x) sum(x^(1:10))) + 1))
 
  
  normal <- mapply(function(x, x_0, d_ii, rho, mu) mvtnorm::dmvnorm(c(x, x_0), c(mu[d_ii], mu[d_ii]), cbind(c(1.53284, rho[d_ii]), c(rho[d_ii], 1.53284)),
                   log = T), x_i, x_menos, d_i, MoreArgs = list(rho = rho, mu = mu))

  return(sum(log(wj)) + sum(normal) + denom)
}
#0.3518519
ls <- label.switching::label.switching(c("ECR-ITERATIVE-1", "AIC", "DATA-BASED"),
                                       z = do.call(rbind, itera[[6]])[90001:100000,][(1:10000)%%10 == 0,],
                                       constraint = 1,
                                       complete = complete.normal.loglikelihood,
                                       mcmc = mcmc.pars, data = o3[-1], K = 15)

mcmc_per <- permute.mcmc(mcmc.pars, ls$permutations$AIC)
coeficientes <- colMeans(mcmc_per[[1]][,,3])
pesos <- colMeans(mcmc_per[[1]][,,2])
mus <- colMeans(mcmc_per[[1]][,,1])

# sim_x <- datos$promedios
sim_x <- o3
for(i in 2:length(sim_x)){
  sim_x[i] <- rmodelo(o3[i-1], tau = 0.652384 ,
                      rho = coeficientes, mus = mus, wj = pesos)
}
sim_x <- x*34.02042
for(i in 2:length(sim_x)){
  sim_x[i] <- rnorm(1,sim_x[i-1],1)
}

# for(i in 2:length(x_n)){
#   x_n[i] <- rmodelo(x_n[i-1], tau = 0.8,
#                       rho = rep(0.99,3), mus = c(-1,0,3), wj = c(0.1,0.4,0.5))
# }



