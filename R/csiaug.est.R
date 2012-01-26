

###Note: Incidence estimation code is at the top;
###Example using simulated Data is at the bottom.

Gupt  =  function(x, k, g)
{
    return(pweibull(x, k, g, FALSE))
}

muFunc = function(k, g)
{
    return(g * gamma(1 + 1 / k))
}

f = function(x, k, g)
{
    return(Gupt(x, k, g) / muFunc(k, g))
}

fkg = function(x, k, g)
{
    return(f(x, k, g))
}

F = function(x, k, g)
{
    return(integrate(function(t){fkg(t, k, g)}, 0, x)$value)
}

f2 = function(x, k, g, mu)
{
    return(Gupt(x, k, g) / mu)
}

fkg2 = function(x, k, g, mu)
{
    f2(x, k, g, mu)
}

F2 = function(l, u, k, g, mu)
{
    integrate(function(t){fkg2(t, k, g, mu)}, l, u)$value

}

L3w = function(params, ai, bi, n1)
{
    k = 20 * exp( - exp(params[1]))
    mu = exp(params[2])
    g = mu / gamma(1 + 1 / k)
    
    first = 0
    for(i in 1:n1)
    {
        first = first + log(max(F2(ai[i], bi[i], k, g, mu), 1e-8, na.rm=TRUE))
    }
    return( - 1 * (first))
}
 
 

##############
csiaug.est <-  function(N1, n1, n0, N3,  ftime, alpha = 0.05) 
{
    ai <- ftime[, 1]
    bi <- ftime[, 2]

    par.est = optim(par = c(log( - log(1 / 20)), log(.5)), fn = L3w, ai = ai, bi = bi, n1 = n1)$par
    
    phat = (N3 + n1) / (N3 + n1 + n0);
    a.hat = 20 * exp( - exp(par.est[1]));
    mu.hat = exp(par.est[2]);
    lambda.hat = n1 / (N1 * mu.hat * phat);
    f.hat = (N3 + n1 + n0) / (N1 + N3 + n1 + n0) * n1 / (N3 + n1) * 1 / mu.hat;
    phi.hat = 1 - (N3 * mu.hat / n1 + mu.hat) * f.hat;
    
     ### calculate individual scores
     ### for those belong to N1 group
    
     U.phi00 = rep(1 / phi.hat, N1);
     U.lambda00 = rep(0, N1);
     U.mu00 = rep(0, N1);
     U.a00 = rep(0, N1);
     U.p00 = rep(0, N1);
    
    ### for those belong to N3 group
    
     denom = 1 - phi.hat - phi.hat * lambda.hat * mu.hat
     U.phi11 = rep( - (1 + lambda.hat * mu.hat) / denom, N3);
     U.lambda11 = rep( - phi.hat * mu.hat / denom, N3);
     U.mu11 = rep( - phi.hat * lambda.hat / denom,  N3);
     U.a11 = rep(0, N3);
     U.p11 = rep(1 / phat, N3);
    
    
    ### for those belong to n1 group
    
     U.phi10 = rep(1 / phi.hat, n1);
     U.lambda10 = rep(1 / lambda.hat,  n1);
     U.mu10 = rep(0, n1);
     U.a10 = rep(0, n1);
    
     tempA = (gamma(1 + 1 / a.hat) / mu.hat)^a.hat;
     tempB = lgamma(1 + 1 / a.hat) - log(mu.hat) - digamma(1 + 1 / a.hat) / a.hat
    
     intA = rep(0, n1);
     intB = rep(0, n1);
     intC = rep(0, n1);
    
     i <- 1
     while (i <= n1) {
         temp.samp = runif(10000, ai[i], bi[i]);
         intA[i] = mean(exp( - temp.samp^a.hat * tempA),  na.rm = TRUE);
         intB[i] = mean(temp.samp^a.hat * exp( - temp.samp^a.hat * tempA),  na.rm = TRUE);
         intC[i] = mean(log(temp.samp) * temp.samp^a.hat * exp( - temp.samp^a.hat * tempA),  na.rm = TRUE);
         U.mu10[i] = (a.hat / mu.hat) * tempA * intB[i] / intA[i];
         U.a10[i] =  - (tempA * intC[i] + tempA * tempB * intB[i]) / intA[i];
         i <- i + 1
    }
    
    U.p10 = rep(1 / phat, n1);
    
    
    if(n0 != 0) {
    
        ### for those belong to n0 group; only calculate this when n0 is not 0;
       
        U.phi.n0 = rep( - 1 / (1 - phi.hat), n0);
        U.lambda.n0 =  rep(0, n0);
        U.mu.n0 = rep(0, n0);
        U.a.n0 = rep(0, n0);
        U.p.n0 = rep( - 1 / (1 - phat), n0);
       
       
        U.phi.old = c(U.phi00, U.phi10, U.phi11, U.phi.n0);
        U.lambda.old = c(U.lambda00, U.lambda10, U.lambda11, U.lambda.n0);
        U.mu.old = c(U.mu00, U.mu10, U.mu11, U.mu.n0);
        U.a.old = c(U.a00, U.a10, U.a11, U.a.n0);
        U.p.old = c(U.p00, U.p10, U.p11, U.p.n0);
       
        U.phi = U.phi.old - mean(U.phi.old,  na.rm = TRUE);
        U.lambda = U.lambda.old - mean(U.lambda.old,  na.rm = TRUE);
        U.mu = U.mu.old - mean(U.mu.old,  na.rm = TRUE);
        U.a = U.a.old - mean(U.a.old,  na.rm = TRUE);
        U.p = U.p.old - mean(U.p.old,  na.rm = TRUE);
       
        Lam = cbind(U.lambda, U.phi, U.mu, U.p);
        MM = t(Lam) %*% Lam - (t(Lam) %*% U.a) %*%(t(U.a) %*%Lam) / (sum(U.a * U.a,  na.rm = TRUE));
       
        cov.est = solve(MM);
        se.est = sqrt(diag(cov.est));
    }
    
    else {
      
        U.phi.old = c(U.phi00, U.phi10, U.phi11);
        U.lambda.old = c(U.lambda00, U.lambda10, U.lambda11);
        U.mu.old = c(U.mu00, U.mu10, U.mu11);
        U.a.old = c(U.a00, U.a10, U.a11);
     
        U.phi = U.phi.old - mean(U.phi.old,  na.rm = TRUE);
        U.lambda = U.lambda.old - mean(U.lambda.old,  na.rm = TRUE);
        U.mu = U.mu.old - mean(U.mu.old,  na.rm = TRUE);
        U.a = U.a.old - mean(U.a.old,  na.rm = TRUE);
     
        Lam = cbind(U.lambda, U.phi, U.mu);
        MM = t(Lam) %*% Lam - (t(Lam) %*% U.a) %*% (t(U.a) %*% Lam) / (sum(U.a * U.a,  na.rm = TRUE));
    
        cov.est = solve(MM);
        se.est = sqrt(diag(cov.est));
    }
    
    parm.est = c(lambda.hat, phi.hat, mu.hat, phat)
    if(length(se.est) == 3)
    {
        se.est[4] = NA
    }
    parms = rbind(parm.est, se.est)
    colnames(parms) = c("lambda", "phi", "mu", "p")
    za = qnorm(1 - alpha / 2)
    lambda.ci = exp(log(lambda.hat) + se.est[1] / lambda.hat * za * c( - 1, 1))
    
    return(list(parameters = round(parms, 5), incidence.ci = round(lambda.ci, 5)))
}


