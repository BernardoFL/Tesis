###########################
#### Full conditionals ####
###########################

#### w_j

vj <- function(j, b, d, z){
  #b es el hiperparámetro
  zi <- unlist(z)
  nj <- length(which(dj == j))
  Nj <- length(which(zi == j))
  nj_p <- length(which(dj > j))
  Nj_p <- length(which(zi > j))
  return(rbeta(1, 1 + nj + Nj, b + nj_p + Nj_p))
}

wj <- function(j, vj, zi, di, b, epsilon){
  zi <- unlist(zi)
  J1 <- max(zi, di)
  inv <- cumprod(1-vj)
  men <- inv[which(inv < epsilon)]
  J2 <- which.min(men)
  J <- max(J1, J2)
  vjs <- sapply(1:J, vj, b = b, d = di, z = zi)
  wjs <- vector('double', J)
  wjs[1] <- vjs[1]
  wjs[2:J] <- vjs[2:j]*cumprod(1 - vjs[1:(J-1)])
  return(wjs)
}

#### \rho_j

rho_j <- function(x, old_rho, R, pi, d_i, sigma, w_j){
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
      return( pi(r)*(1-r^2)^(-length(which(dj == j)) / 2) * exp(-tau/2*suma) )
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

#### d_i


d_i <- function(x, old_d, wj, ksi, mu, sigma, rho){
  vj <- sapply(old_d, function(x) runif(1, min = 0, max = exp(-ksi*x)) )
  Js <- sapply(vj, function(x) floor( -(1/ksi)*log(x) ) )
  sop <- lapply(Js, function(x) return(1:x))
  dens <- function(x, y, mu, mu.y, rho){ function(j){exp(ksi*j)*wj[j]*dnorm(x, mu, sigma)*dnorm(y, mu.y+rho*(x-mu.y), sigma) }}
  app.dens <- vector('list', length(x)-1)
  for(i in 2:length(x)){
    app.dens[[i-1]] <- dens(x[i], x[i-1], mu[i], mu[i-1], rho[i]) 
    #dens(x[i], x[i-1], mu, mu.y, rho) #Necesito checar los parámetros
  }
  app.dens[[1]] <- dens(x[2], x[1], mu[2], mu[1], rho[2])
  probs <- vector('list', length(sop))
  for(i in 1:length(sop)){
    probs[[i]] <- sapply(sop[[i]], app.dens[[i]])
  }
  probs <- lapply(probs, function(x) x/sum(x))
  d.sample <- function(pr) sample(1:length(pr), 1, prob = pr)
  new.d <- lapply(probs, d.sample)
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
  for(i in 2:length(x)){
    app.dens[[i-1]] <- dens(x[i-1], mu[i]) 
  }
  app.dens[[1]] <- dens(x[1], mu)#Ni idea por qué pero sólo así jala
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

mu_j <- function(j, x, old_mu, rho, d_i, z_il, t, sigma){
  t_j <- function(j){ t + 2* length(which(d_i == j)) /(1+rho[j])  }
  m_j <- function(j, t.j){
    indices <- which(d_i == j)
    suma <- sum(x[indices] + x[indices-1])/(1+rho[j])
    return(  1/t.j* (m*t + (1/sigma) * suma ) )
  }
  num.prob <- vector(length(x)-1, 'double')
  for(i in 1:length(x)){
    num.prob[i] <- (x[i-1] -  old_mu[ z_il[[i]] ])^2
  } 
  num.prob <-  1 - exp( -num.prob/(2*sigma) )
  num.prob <- prod(num.prob)

  
  mh.step <- function(j){
    t.j <- t_j(j)
    mu.prop <- rnorm(1, m_j(j, t.j), 1/t.j)
    indices <- unlist(lapply(z_il, function(x) length(which(x == j)) ))
    denom.uno <-  (1 - exp( (x[which(indices > 0) - 1] - mu.prop)^2 / (2*sigma) ))^indices[which(indices > 0)] 
    indices.dos <- unlist(lapply(z_il, function(x) length(which(x != j)) ))
    denom.dos <- (1 - exp( (x[which(indices.dos > 0) - 1] - mu.prop)^2 / (2*sigma) ))^indices[which(indices.dos > 0)] 
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

tau <- function(x, mu, d_i, old_tau, rho, c, a, z_il){
  require(truncdist)
  #c
  sigma_r <- function(r)  return(solve( cbind(c(1, r), c(r, 1)) ))
  mu_i <- function(i) return( t( c(x[i+1] - mu[d_i[i]], x[i] - mu[d_i[i]] ) ) )
  mus <- sapply(1:length(x) - 1, mu_i)
  s_r <- sapply(rho, sigma_r)
  aux <- function(x) t(x) %*% s_r %*% x
  suma <- sapply(mus, aux)
  suma <- sum(suma)
  
  c.hat <- c + suma/2
  a.hat <- a + length(d_i)/2 
  
  aux <- vector(length(z_il), 'list')
  for(i in 1:length(z_il)){
    aux[[i]] <- x[i] - mu[ z_il[[i]] ]
  }
  aux <- unlist(aux)
  aux <- exp( -old_tau * aux^2/2 )
  u_i <- sapply(aux, function(x) runif(1, x, 1))
  aux.dos <- vector(length(z_il), 'list')
  for(i in 1:length(z_il)){
    aux.dos[[i]] <- x[i+1] - mu[  z_il[[i]]  ]
  }
  aux.dos <- unlist(aux.dos)
  te <- max(-2*log(1 - u_i)/aux.dos^2)
  new.tau <- rtrunc(1, 'gamma', te, Inf)
  return(1/new.tau)
}


### k_i

k_i <- function(old_k, x, wj, p, z_il){
  z.new <- sample(1:length(wj), length(old_k), T, prob = wj)
  comp <- vector(length(x), 'double')
  for(i in 1:length(x)){
    comp[i] <- x[i+1] - mu[z.new[i]]
  }
  prob <- p/(1-p) * (1 - exp(-1/(2*sigma) * comp^2) )
  prob <- sapply(prob, function(x) min(1, x))
  ratios <- runif(length(prob))
  
  comp.des <- vector(length(x), 'double')
  for(i in 1:length(x)){
    comp.des[i] <- x[i+1] - mu[  z_il[[i]][old_k[i]]  ]
  }
  prob.des <- (1-p)/p * 1/(1 - exp(-1/(2*sigma) * comp.des^2) )
  prob.des <- sapply(prob.des, function(x) min(1, x))
  
  no.acepto <- sapply(1:length(prob), function(i) ifelse(ratios[i] > prob[i], TRUE, FALSE))
  new.k <- vector(length(old_k), 'numeric')
  new.k <- sapply(no.acepto, function(i) ifelse(no.acepto[i], new.k[i] <- old_k[i],
                                                new.k[i] <- old_k[i] + 1) )
  descend <- function(i, no.acepto){
    if(no.acepto){
      r <- runif(1)
      ifelse(r < prob.des[i], new.k[i] <- old.k[i] - 1, new.k[i] <- old.k[i])
    }
  }
  return(new.k)
}


##########


simulacion <- function(n){
  dens <- function(x) 0.2*dnorm(x, -5, 1)+0.4*dnorm(x, 0, 3)+0.4*dnorm(x, 7) 
  dens <- Vectorize(dens, 'x')
  curve(dens, xlim=c(-15,15))
  
  componentes <- sample(1:3,prob=c(0.2, 0.4, 0.4),size=n,replace=TRUE)
  mus <- c(-5, 0, 7)
  sds <- sqrt(c(1, 3, 1))
  
  samples <- rnorm(n, mean=mus[componentes], sd=sds[componentes])
  return(samples)
}

post.inf <- function(datos, nsim){
  simul <- vector(nsim, 'list')
  init 
  for(i in 1:nsim){
    simul[[i]] <- list()
  }
}