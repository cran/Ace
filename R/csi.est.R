###Cross-Sectional Incidence####
###Parameters (mu,p) assumed to be known precisely by default###


csi.est = function(N1, N2, N3, mu = 0.5, p = 1, se.mu = 0, se.p = 0, alpha = 0.05)
{
    inc.est = (p * N2 - (1 - p) * N3) / (p * N1 * mu)
    a = N2 / (N1^2 * mu^2 * p)
    b = (se.p^2 * N3^2) / (N1^2 * mu^2 * p^4)	
    c = (se.mu^2 * inc.est^2) / mu^2
    var.est = a + b + c
    za = qnorm(1 - alpha / 2)
    ci = inc.est * exp(c(-1,1) * za * sqrt(var.est) / inc.est)
    res <- list(inc.est = inc.est, se.est = sqrt(var.est), incidence.ci = ci)
    return(res)
}
    
