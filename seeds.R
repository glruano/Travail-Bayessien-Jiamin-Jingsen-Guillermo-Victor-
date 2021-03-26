library(invgamma)

pi <- function(x1, x2, alpha0, alpha1, alpha2, alpha12, b){
  return(exp(alpha0*rep(1,length(x1)) + alpha1*x1 + alpha2*x2 + alpha12*x1*x2 + b)/(exp(alpha0*rep(1,length(x1)) + alpha1*x1 + alpha2*x2 + alpha12*x1*x2 + b)+1))
}


seeds <- function(n, N, x1, x2, nchain = 10^4, prop_sd = c(0.1, 0.2, 0.15, 0.3, 0.3)){
  I <- length(n)
  chain <- matrix(NA, nchain + 1, 5+I)
  colnames(chain) <- c("alpha0", "alpha1", "alpha2", "alpha12", paste("b", 1:I, sep =""), "sigma2")
  #initialization
  chain[1,] <- rep(0, 5+I)
  acc_rates <- rep(0, 4+I)
  
  #alpha0, alpha1, alpha2, alpha12, bi, sigma2
  
  for (iter in 1:nchain){
    current <- chain[iter,]
    
    ##Mise à jour sigma2
    update_shape <- 10^(-3) + I/2
    update_scale <- 10^(-3) + sum(current[5:(4+I)]^2)/2
    sigma2 <- rinvgamma(1, shape=update_shape, scale=update_scale)
    current[5+I] <- sigma2
 
    ## Mise a jour de alpha0
    prop <- current
    prop[1] <- rnorm(1, current[1], prop_sd[1])
    top <- -0.5*prop[1]^2/10^6 + sum(n*log(pi(x1, x2, prop[1], current[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, prop[1], current[2], current[3], current[4], current[5:(4+I)])))
    bottom <- -0.5*current[1]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)])))
    acc_prob <- exp(top-bottom)
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[1] <- acc_rates[1] + 1
    }
    
    ## Mise a jour de alpha1
    prop <- current
    prop[2] <- rnorm(1, current[2], prop_sd[2])
    
    top <- -0.5*prop[2]^2/10^6 + sum(n*log(pi(x1, x2, current[1], prop[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], prop[2], current[3], current[4], current[5:(4+I)])))
    bottom <- -0.5*current[2]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)])))
    acc_prob <- exp(top-bottom)
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[2] <- acc_rates[2] + 1
    }
    
    
    ## Mise a jour de alpha2
    prop <- current
    prop[3] <- rnorm(1, current[3], prop_sd[3])
    
    top <- -0.5*prop[3]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], prop[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], prop[3], current[4], current[5:(4+I)])))
    bottom <- -0.5*current[3]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)])))
    acc_prob <- exp(top-bottom)
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[3] <- acc_rates[3] + 1
    }
    
    ## Mise a jour de alpha12
    prop <- current
    prop[4] <- rnorm(1, current[4], prop_sd[4])
    
    top <- -0.5*prop[4]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], current[3], prop[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], current[3], prop[4], current[5:(4+I)])))
    bottom <- -0.5*current[4]^2/10^6 + sum(n*log(pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)]))+(N-n)*log(rep(1,I)-pi(x1, x2, current[1], current[2], current[3], current[4], current[5:(4+I)])))
    acc_prob <- exp(top-bottom)
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[4] <- acc_rates[4] + 1
    }
    
    ## Mise a jour de b
    
    for (i in 1:I){
      prop <- current
      prop[4+i] <- rnorm(1, current[4+i], prop_sd[5])
      top <- -0.5*prop[4+i]^2/current[5+I]^2 + n[i]*log(pi(x1[i], x2[i], current[1], current[2], current[3], current[4], prop[4+i]))+(N[i]-n[i])*log(1-pi(x1[i], x2[i], current[1], current[2], current[3], current[4], prop[4+i]))
      bottom <- -0.5*current[4+i]^2/current[5+I]^2 + n[i]*log(pi(x1[i], x2[i], current[1], current[2], current[3], current[4], current[4+i]))+(N[i]-n[i])*log(1-pi(x1[i], x2[i], current[1], current[2], current[3], current[4], current[4+i]))
      acc_prob <- exp(top-bottom)
      
      if (runif(1) < acc_prob){
        current <- prop
        acc_rates[4+i] <- acc_rates[4+i] + 1
      }
    }
    
    
    ## Sauvegardons le nouvel etat
    chain[iter+1,] <- current
  }
  return(list(chain = chain, acc_rates = acc_rates / nchain))
}

## Application
I <- 21
n <- c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 
    3, 22, 15, 32, 3)
N <- c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 
    4, 12, 41, 30, 51, 7)
x1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1)
x2 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 
    1)

out <- seeds(n, N, x1, x2)
out$chain <- out$chain[-(1:1000),]

par(mfrow = c(5, 2), mar = c(4, 5, 0.5, 0.5))
for (j in 1:4){
  plot(out$chain[,j], type = "l", main = "")
  plot(density(out$chain[,j]), type = "l", main = "")
}

plot(out$chain[,26], type = "l", main = "")
plot(density(out$chain[,26]), type = "l", main = "")


moy_alpha0 <- mean(out$chain[,1])
moy_alpha1 <- mean(out$chain[,2])
moy_alpha2 <- mean(out$chain[,3])
moy_alpha12 <- mean(out$chain[,4])
moy_sigma <- sqrt(mean(out$chain[,26]))
