###########################
#### Full conditionals ####
###########################

#### w_j

v_j <- function(j, b, d, z){
  #b es el hiperparámetro
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
  vjs <- sapply(1:(J1*5), v_j, b = b, d = di, z = zi)
  inv <- cumprod(1-vjs)
  men <- inv[which(inv < epsilon)]
  J2 <- which.min(men)
  J <- max(J1, J2)
  vjs <- vjs[1:J]
  wjs <- vector('double', J)
  wjs[1] <- vjs[1]
  wjs[2:J] <- vjs[2:J]*cumprod(1 - vjs[1:(J-1)])
  if(any(is.na(wjs)))
    message('Hay NAs en los pesos')
  return(wjs)
}

#### \rho_j

rho_j <- function(x, old_rho, R, d_i, mu, sigma, w_j, same_rho){
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
    recover()
    probs <- non_normal / sum(non_normal) 
    probs <- asNumeric(probs)
    return(sample(R, 1, prob = probs))
  } else{
    tau <- 1/sigma
    sigma_r <- function(r)  return(solve( cbind(c(1, r), c(r, 1)) ))
    mu_i <- function(i) return( t( c(x[i+1] - mu[d_i[i]], x[i] - mu[d_i[i]] ) ) )
    dens <- function(j){
      function(r){
        s_r <- sigma_r(r)
        sumable <- which(d_i == j)
        mus <- sapply(sumable, mu_i)
        aux <- function(x) t(x) %*% s_r %*% x
        suma <- sapply(mus, aux)
        suma <- sum(suma)
        return( (1-r^2)^(-length(which(dj == j)) / 2) * exp(-tau/2*suma) )
      }
    }
    dens.app <- lapply(1:length(w_j), dens)
    probs <- matrix(nrow = length(R), ncol = length(w_j))
    for(i in 1:length(dens.app)){
      probs[,i] <- sapply(R, dens[[i]])
    }
    probs <- apply(probs, 2, function(x) x/sum(x))
    new_rho <- apply(probs, 2, function(x) sample(R, 1, prob = x))
    return(new_rho)
  }
}

#### d_i


d_i <- function(x, old_d, wj, ksi, mu, sigma, rho, same_rho, i){
  sigma <- sqrt(sigma)
  vj <- sapply(old_d, function(x) runif(1, min = 0, max = exp(-ksi*x)) )
  Js <- sapply(vj, function(x) floor( -(1/ksi)*log(x) ) )
  sop <- lapply(Js, function(x) return(1:x))
  dens <- function(x, y, mu, mu.y, rho){ function(j){exp(ksi*j)*wj[j]*dnorm(mpfr(x,100), mu, mpfr(sigma,100))*dnorm(mpfr(y,100), mu.y+rho*(x-mu.y), mpfr(sigma,100)) }}
  app.dens <- vector('list', length(x)-1)
  for(i in 2:(length(x) - 1)){#checar este ciclo
    if(same_rho) #n+1
      app.dens[[i]] <- dens(x[i+1], x[i], mu[i], mu[i-1], rho) 
    else
      app.dens[[i]] <- dens(x[i+1], x[i], mu[i], mu[i-1], rho[i])       
    #dens(x[i], x[i-1], mu, mu.y, rho) #Necesito checar los parámetros
  }
  if(same_rho)
    app.dens[[1]] <- dens(x[2], x[1], mu[2], mu[1], rho)    
  else
    app.dens[[1]] <- dens(x[2], x[1], mu[2], mu[1], rho[2])
  probs <- vector('list', length(sop))
  for(i in 1:length(sop)){
    probs[[i]] <- do.call(app.dens[[i]], list(sop[[i]]))
  }
  probs <- lapply(probs, function(x) x/sum(x))
  probs <- lapply(probs, asNumeric)
  d.sample <- function(pr) sample(1:length(pr), 1, prob = pr)
  new.d <- lapply(probs, function(x) sample(1:length(x), 1, prob = x))
  new.d <- unlist(new.d)
  return(new.d)
}


#### z_li

z_li <- function(x, old_z, wj, phi, mu, sigma, ki){
  runifo <- Vectorize(runif, 'max')
  vj <- lapply(old_z, function(x) runifo(1, min = 0, max = exp(-phi*x) ) )
  Js <- lapply(vj, function(x) sapply(x, function(y) floor( -(1/phi)*log(y) ) ) )
  sop <- lapply(Js, function(x) lapply(x, function(y) return(1:y))  )
  dens <- function(y, mu){ function(j){exp(phi*j)*wj[j]*(1-exp(-(y-mu)^2/(2*sigma))) }}
  app.dens <- vector('list', length(x)-1)
  for(i in 2:(length(x) - 1)){
    app.dens[[i]] <- dens(x[i], mu[i]) 
  }
  app.dens[[1]] <- dens(x[1], mu[1])#Ni idea por qué pero sólo así jala
  probs <- vector('list', length(sop))
  for(i in 1:length(sop)){
    probs[[i]] <- lapply(sop[[i]], function(x) sapply(x, app.dens[[i]]) )
  }
  probs <- lapply(probs, function(x) lapply(x, function(y) y/sum(y)) )
  z.sample <- function(pr) sample(1:length(pr), 1, prob = pr)
  new.z <- lapply(probs, function(x) lapply(x, z.sample) )
  new.z <- lapply(new.z, unlist)
  return(new.z)
}


#### mu_j

mu_j <- function(j, x, old_mu, rho, d_i, z_il, t, m, sigma, same_rho){
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
  
  num.prob <- vector('list', length(x)-1)
  for(i in 1:(length(x) - 1)){
    num.prob[[i]] <- 1 - exp( -(x[i] -  old_mu[ z_il[[i]] ])^2/(2*sigma) )
  } 
  num.prob <-  unlist(num.prob)
  num.prob <- prod(num.prob)

  mh.step <- function(j){
    t.j <- t_j(j)
    mu.prop <- rnorm(1, m_j(j, t.j), 1/t.j)
    indices <- unlist(lapply(z_il, function(x) length(which(x == j)) ))
    denom.uno <-  (1 - exp( (x[which(indices > 0) ] - mu.prop)^2 / (2*sigma) ))^indices[which(indices > 0)] 
    denom.uno <- prod(denom.uno)
    indices.dos <- unlist(lapply(z_il, function(x) length(which(x != j)) ))
    denom.dos <- (1 - exp( (x[which(indices.dos > 0)] - mu.prop)^2 / (2*sigma) ))^indices[which(indices.dos > 0)] 
    denom.dos <- prod(denom.dos)
    prob <- num.prob / (denom.uno * denom.dos)
    prob <- min(1, prob)
    ratio <- runif(1)
    ifelse(ratio < prob, return(mu.prop), return(old_mu[j]))
  }
  ## Checar este pedo
  new.mu <- sapply(1:length(old_mu), mh.step)
  return(new.mu)
}


#### tau

tau <- function(x, mu, d_i, old_tau, rho, c, a, z_il, same_rho){
  sigma_r <- function(r)  solve( cbind(c(1, r), c(r, 1)) )
  mu_i <- function(i) return( t( c(x[i+1] - mu[d_i[i]], x[i] - mu[d_i[i]] ) ) )
  mus <- sapply(1:(length(x) - 1), mu_i)
  if(same_rho){
    s_r <- sigma_r(rho)
    aux <- function(i) t(mus[,i]) %*% s_r %*% mus[,i]
  } else{
    s_r <- sapply(rho, sigma_r)
    aux <- function(i) t(mus[,i]) %*% s_r %*% mus[,i]    
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
  aux <- exp( -aux^2/(2*old_tau) )
  u_i <- sapply(aux, function(x) runif(1, x, 1))
  aux.dos <- vector('list', length(z_il))
  for(i in 1:length(z_il)){
    aux.dos[[i]] <- x[i+1] - mu[  z_il[[i]]  ]
  }
  aux.dos <- unlist(aux.dos)
  te <- mpfr(max(-2*log(1 - u_i)/aux.dos^2), prec = 100)
  new.tau <- rtruncgamma(a.hat, c.hat, old_tau, te)
  if(1/new.tau == 0)
    message('Underflow en la varianza')
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
  recover()
  for(i in 1:(length(x)-1)){
    comp.des[i] <- x[i+1] - mu[  z_il[[i]][old_k[i]]  ]
  }
  prob.des <- (1-p)/p * 1/(1 - exp(-1/(2*sigma) * comp.des^2) )
  prob.des <- sapply(prob.des, function(x) min(1, x))
  
  no.acepto <- sapply(1:length(prob), function(i) ratios[i] > prob[i])
  new.k <- numeric(length(old_k))
  for(i in 1:length(no.acepto)){
    if(no.acepto[i]){
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
  eval(parse(text=paste("z[[i]] <- list(", paste(z_il, collapse=','), ')')),
             envir = parent.frame())# :O
  return(new.k)
}

rtruncgamma <- function(shape, rate, prev.tau, trunc){
  trunc <- mpfr(trunc, prec = 100); rate <- mpfr(rate, prec=100)
  y <- mpfr(runif(1, 0, exp(-prev.tau)), prec = 100)
  x <- mpfr(runif(1), prec = 100)
  (((x*( (-log(y))^shape )) + ((1-x)*trunc^shape))^(1/shape))/rate
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
mcmc <- function(datos, nsim, ksi, phi, R, same_rho, b, c, a, p, t, m, epsilon, init){
  suppressPackageStartupMessages(require(Rmpfr))
  w <- vector('list', nsim); w[[1]] <- init[['w']]
  d <- vector('list', nsim); d[[1]] <- init[['d']]
  z <- vector('list', nsim); z[[1]] <- init[['z']]
  rho <- vector('list', nsim); rho[[1]] <- init[['rho']]
  mu <- vector('list', nsim); mu[[1]] <- init[['mu']]
  sigma <- vector('list', nsim); sigma[[1]] <- init[['sigma']]
  k <- vector('list', nsim); k[[1]] <- init[['k']]
    #list(R, mu, dj, z_il, sigma)
  for(i in 2:nsim){
  #checar si el razonamiento es correcto
    d[[i]] <- d_i(x = datos, old_d = d[[i-1]], 
            wj = w[[i-1]], ksi=ksi, mu = mu[[i-1]],
            sigma = sigma[[i-1]], rho = rho[[i-1]], same_rho = same_rho, i=i)
    z[[i]] <- z_li(x = datos, old_z = z[[i-1]], wj = w[[i-1]], 
               phi=phi, mu = mu[[i-1]], sigma = sigma[[i-1]], ki = k[[i]])
    w[[i]] <- wj(z[[i]], d[[i]], b, epsilon)
    rho[[i]] <- rho_j(x = datos, old_rho = rho[[i-1]], R,
                d_i = d[[i]], mu = mu[[i-1]],
                sigma = sigma[[i-1]], w_j = w[[i]], same_rho = same_rho)
    mu[[i]] <- mu_j(j = i, x = datos, old_mu = mu[[i-1]],
              rho = rho[[i]], d_i = d[[i]], m = m,
              z_il = z[[i]], t, sigma = sigma[[i-1]], same_rho=same_rho)
    sigma[[i]] <- tau(x = datos, mu = mu[[i]], d_i = d[[i]],
              old_tau = 1/sigma[[i-1]], rho = rho[[i]],
              c, a, z_il = z[[i]], same_rho = same_rho)
    k[[i]] = k_i(old_k = k[[i-1]], x = datos,
            wj = w[[i]], p, z_il = z[[i]], mu = mu[[i]], sigma = sigma[[i]])
  }
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
             rho = 0.5,
             sigma = 3,
             mu = sample(c(-1,0.5,4), 800, replace = T))  #ver si estas madres son convexas
#nota: checar x_0
itera <- mcmc(x_n, 100, ksi = 1, phi = 1,
              R = seq(0.001, 0.999, by = 0.001), same_rho = T,
              a = 1, c = 0.1, p = 0.5, m = mean(x_n), b = 0.1,
              t = 1/sd(x_n), epsilon = 0.1, init = init)

