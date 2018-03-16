###########################
#### Full conditionals ####
###########################

#### w_j

v_j <- function(j, b, d, z){
  zi <- unlist(z)
  nj <- length(which(d == j))
  Nj <- length(which(zi == j))
  nj_p <- length(which(d > j))
  Nj_p <- length(which(zi > j))
  return(rbeta(1, 1 + nj + Nj, b + nj_p + Nj_p))
}

wj <- function(zi, di, b, epsilon){
  zi <- unlist(zi)
  J1 <- max(zi, di)
  vjs <- sapply(1:(J1*40), v_j, b = b, d = di, z = zi)
  inv <- cumprod(1-vjs)
  men <- which(inv < epsilon)
  J2 <- which.min(men)
  J <- max(J1, J2)
  vjs <- vjs[1:J]
  wjs <- vector('double', length(vjs))
  wjs[1] <- vjs[1]
  wjs[-1] <- vjs[-1]*cumprod(1 - vjs[-length(vjs)])
  if(any(is.na(wjs)))
    message('Hay NAs en los pesos')
  return(wjs)
}

#### \rho_j 

rho_j <- function(x, old_rho, R, d_i, mu, sigma, w_j, same_rho, error){
  if(same_rho){
    non_normal <- mpfrArray(0, prec = 100, dim = c(1,length(R)))
    dens <- function(r){
      x_prev <- c(x[-c(1, length(x))],0)
      expo <- -sum((x[-1] - mu[d_i]*(1-r) - r*x_prev)^2)/(2*sigma*(1-r^2))
      (mpfr(1-r^2, 100)^(-(length(x) - 1)/2))*exp(mpfr(expo, 100))# n-1
    }
    for(i in 1:length(non_normal)){
      non_normal[i] <- dens(R[i])
    }

    probs <- non_normal / sum(non_normal) 
    probs <- asNumeric(probs)
    return(sample(R, 1, prob = probs))
  } else{
    tau <- 1/sigma
    sigma_r <- function(r)  return(solve( cbind(c(1, r), c(r, 1)) ))
    mu_i <- function(i) return( t( c(x[i+1] - mu[d_i[i]], x[i] - mu[d_i[i]] ) ) )
    dens <- function(j, r){
        s_r <- sigma_r(r)
        sumable <- which(d_i == j)
        if(identical(sumable, integer(0))) #si el cluster estaa vacio
          return(0)
        mus <- sapply(sumable, mu_i)
        aux <- function(x) t(x) %*% s_r %*% x
        suma <- apply(mus, 2, aux)
        suma <- sum(suma)
        
        ((-length(which(d_i == j)) / 2)*log(1-r^2))  - tau*suma/2 
    }
  
    probs <- matrix(nrow = length(R), ncol = length(w_j))
    for(i in 1:ncol(probs)){
      probs[,i] <- vapply(R, dens, 0, j = i)
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
    
    apply(probs, 2, function(x) sample(R, 1, prob = x))
  }
}

#### d_i


d_i <- function(x, old_d, wj, mu, sigma, rho, same_rho, i, error){
  sigma <- sqrt(sigma); n <- length(x) - 1
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
    sorted <- sorted - sorted[length(sorted)]
    sorted <- ifelse(sorted >= log(error) - log(length(sorted)), exp(sorted), 0)
    norms <- sorted/sum(sorted)
    norms[indices]
  })
  apply(probs, 2, function(x) sample(1:length(wj), 1, prob = x))
}


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

mu_j <- function(j, x, old_mu, rho, d_i, z_il, t, m, sigma, same_rho, wj){
  if(same_rho){
    t_j <- function(j){ t + 2* length(which(d_i == j)) /(1+rho)  }    
  } else{
    t_j <- function(j){ t + 2* length(which(d_i == j)) /(1+rho[j])  }
  }  
  m_j <- function(j, t.j){
    indices <- which(d_i == j)
    if(same_rho){
      suma <- sum(x[indices+1] + x[indices])/(1+rho)      
    } else{
      suma <- sum(x[indices+1] + x[indices])/(1+rho[j])
    }
    return(  1/t.j* (m*t + (1/sigma) * suma ) )
  }
  
  old_mu <- old_mu[1:length(wj)]

  mh.step <- function(j){
    num.prob <- vector('list', length(x)-1)
    for(i in 1:(length(x) - 1)){
      num.prob[[i]] <- 1 - exp( -(x[i] -  old_mu[ z_il[[i]] ])^2/(2*sigma) )
    } 
    num.prob <-  unlist(num.prob)
    num.prob <- sum(log(num.prob))
    
    
    
    t.j <- t_j(j)
    mu.prop <- rnorm(1, m_j(j, t.j), sqrt(1/t.j))
    indices <- unlist(lapply(z_il, function(x) length(which(x == j)) ))
    if( identical( which(indices > 0), integer(0) )){
      denom.uno <- 0
      indices <- rep(0, length(indices))
    }
    else{
      exponente <- indices[which(indices > 0)] #es para cuando mas de una vez hace match el z_il 
      if(identical(exponente, numeric(0))) 
        exponente <- 0
      denom.uno <-  (1 - exp( -(x[which(indices > 0) + 1] - mu.prop)^2 / (2*sigma) ))^exponente 
      
      denom.uno <- sum(log(denom.uno))
    }
    

    indices.dos <- unlist(lapply(z_il, function(x) length(which(x != j)) ))
    if( identical( indices.dos, integer(0) ) ){
      denom.dos <- 0 
      indices.dos <- rep(0, length(indices.dos))
    }
    else{
      exponente <- indices[which(indices.dos > 0)] 
      if(identical(exponente, numeric(0))) 
        exponente <- 0
      
      denom.dos <- Map(function(z, equis, poten)  (1 - exp( -(equis - old_mu[unlist(z)])^2 / (2*sigma)) )^poten,
          z_il[which(indices.dos > 0)], x[which(indices.dos > 0)+1], exponente )
      
      denom.dos <- unlist(denom.dos)
      denom.dos <- sum(log(denom.dos))
    } 
  
   
    prob <- exp(num.prob - denom.uno - denom.dos)
    prob <- min(1, prob)
    ratio <- runif(1)
    if(ratio < prob)
      return(mu.prop)
    else
      return(old_mu[j])
  }
  
  for(j in 1:length(wj)){
    old_mu[j] <- mh.step(j)
  }
  return(old_mu)
}


#### tau

tau <- function(x, mu, d_i, old_tau, rho, c, a, z_il, same_rho){
  sigma_r <- function(r)  solve( cbind(c(1, r), c(r, 1)) )
  mu_i <- function(i)  t( c(x[i+1] - mu[d_i[i]], x[i] - mu[d_i[i]] ) ) 
  mus <- sapply(1:(length(x) - 1), mu_i)
  if(same_rho){
    s_r <- sigma_r(rho)
    aux <- function(i) t(mus[,i]) %*% s_r %*% mus[,i]
  } else{
    s_r <- lapply(rho, sigma_r)
    aux <- function(i) t(mus[,i]) %*% s_r[[d_i[i]]] %*% mus[,i]    
  }
  suma <- sapply(1:ncol(mus), aux)
  suma <- sum(suma)
  
  c.hat <- c + suma/2
  a.hat <- a + length(d_i)/2 
  
  aux <- vector('list', length(z_il))
  for(i in 1:length(z_il)){
    aux[[i]] <- x[i] - mu[ z_il[[i]] ]
  }
  aux <- unlist(aux)
  aux <- 1 - exp( -aux^2/(2*old_tau) )
  u_i <- sapply(aux, function(x) runif(1, 0, x))
  aux.dos <- vector('list', length(z_il))
  for(i in 1:length(z_il)){
    aux.dos[[i]] <- x[i+1] - mu[  z_il[[i]]  ]
  }
  aux.dos <- unlist(aux.dos)
  te <- max(-2*log(1 - u_i)/mpfr(aux.dos^2, 100))
  if(te > 10)
    te <- 10
  new.tau <- rtruncgamma(a.hat, c.hat, 1/old_tau, te)
  if(1/new.tau == 0){
    stop('Underflow en la varianza')
  }

  return(asNumeric(1/new.tau))
}


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

rtruncgamma <- function(shape, rate, prev.tau, trunc){
  trunc <- mpfr(trunc, prec = 100); rate <- mpfr(rate, prec=100)
  y <- runif(1)*exp(-rate*mpfr(prev.tau, 100))
  if(y == 0)
    y <- mpfr(1e-320, 100)
  x <- mpfr(runif(1), prec = 100)
  (x * (trunc^shape + (-log(y))^shape )/ (shape*(shape-1))  )^(1/shape)
}
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
mcmc <- function(datos, nsim,  R, same_rho, b, c, a, p, t, m, epsilon, init, error){
  suppressPackageStartupMessages(require(Rmpfr))
  w <- vector('list', nsim); w[[1]] <- init[['w']]
  d <- init[['d']]
  z <- init[['z']]
  rho <- vector('list', nsim); rho[[1]] <- init[['rho']]
  mu <- vector('list', nsim); mu[[1]] <- init[['mu']]
  sigma <- vector('list', nsim); sigma[[1]] <- init[['sigma']]
  k <- vector('list', nsim); k[[1]] <- init[['k']]
 pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    #list(R, mu, dj, z_il, sigma)
  for(i in 2:nsim){
    #print(i)
   if(i %% 100 == 0)
     setTxtProgressBar(pb, i)
    w[[i]] <- wj(z, d, b, epsilon)
    d <- d_i(x = datos, old_d = d, 
            wj = w[[i]], mu = mu[[i-1]],
            sigma = sigma[[i-1]], rho = rho[[i-1]], same_rho = same_rho, i=i, error=error)
    z <- z_li(x = datos, old_z = z, wj = w[[i]], 
              mu = mu[[i-1]], sigma = sigma[[i-1]], ki = k[[i-1]], di = d, error = error)
    rho[[i]] <- rho_j(x = datos, old_rho = rho[[i-1]], R,
                d_i = d, mu = mu[[i-1]],
                sigma = sigma[[i-1]], w_j = w[[i]], same_rho = same_rho, error=error)
    mu[[i]] <- mu_j(j = i, x = datos, old_mu = mu[[i-1]],
              rho = rho[[i]], d_i = d, m = m, wj = w[[i]],
              z_il = z, t, sigma = sigma[[i-1]], same_rho=same_rho)
    sigma[[i]] <- tau(x = datos, mu = mu[[i]], d_i = d,
          c, a, z_il = z, same_rho = same_rho,
old_tau = sigma[[i-1]], rho = rho[[i]])
    k[[i]] = k_i(old_k = k[[i-1]], x = datos,
            wj = w[[i]], p, z_il = z, mu = mu[[i]], sigma = sigma[[i]])
  }
 close(pb)
  return(list(w, rho, mu, sigma))
}

post.inf <- function(mcmc, burn_in){
  w.hat <- colMeans(do.call(cbind, mcmc[['w']])[,-(1:burnin)])
  m.hat <- colMeans(do.call(cbind, mcmc[['mu']])[,-(1:burnin)])
  sigma.hat <- colMeans(do.call(cbind, mcmc[['sigma']])[,-(1:burnin)])
  rho.hat <- colMeans(do.call(cbind, mcmc[['rho']])[,-(1:burnin)])
  return(list(w.hat, m.hat, sigma.hat, rho.hat))
}

paraZ <- sample(1:5, length(x_n) - 1, TRUE)
w <- rbeta(800,1,1)
init <- list(d = sample(1:4, length(x_n)-1, T),
             z = lapply(paraZ, function(n) sample(1:n, n) ),
             w = w/sum(w),
             k = paraZ,
             rho = rep(0.5, 800),
             sigma = 10,
             mu = rep(0.01, 800))  
#nota: checar x_0
itera <- mcmc(x_n, 10000, ksi = 2, phi = 2,
              R = seq(0.001, 0.999, by = 0.01), same_rho = F,
              a = 0.5, c = 1, p = 0.4, m = 0, b = 1,
              t = 1/4, epsilon = 0.00000001, init = init, error=1e-50) 



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



