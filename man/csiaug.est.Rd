\name{csiaug.est}
\alias{csiaug.est}
\title{
Augmented Cross-sectional Estimation of incidence rate
}
\description{
Estimating incidence rate based on augmented cross-sectional studies.
}
\usage{
csiaug.est(N1,  n1,  n0,  N3,  ftime, alpha = 0.05)
}
% -  maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{
Number of subjects in cross-sectional sample testing negative on sensitive test (e.g., ELISA); (Uninfected).
}
  \item{n1}{
Number of subjects in cross-sectional sample testing positive on sensitive test and negative on less-sensitive test (e.g., BED)  AND will eventually test positive on less-sensitive test; (Recently Infected).
}
  \item{n0}{
Number of subjects in cross-sectional sample testing positive on sensitive test and negative on less-sensitive test  AND will continue to test negative on less-sensitive test indefinitely; (Long-term Non-progressors).
}
  \item{N3}{
Number of subjects in cross-sectional sample testing positive on sensitive and less-sensitive; (Non-recently Infected).
}
  \item{ftime}{
\eqn{n_1 \times 2} matrix giving interval-censored observation of forward recurrence time 
(time from cross-sectional sample to first positive result on less-sensitive test or end of follow-up, whichever occurs first).
 Times should be expressed on the same scale for which incidence rate is desired (e.g. years, in order to obtain annual incidence rate).
}
  \item{alpha}{
specifys the coverage level (\code{100*(1-alpha)\%} ) of resulting confidence interval for incidence rate.
}
}
\details{
If no upper bound for \code{ftime} is known, do not use Inf. A finite upper bound is required. A value of 20 years is generally sufficiently large if no information is known.
}
\value{
\item{parameters }{A \eqn{2 \times 4} matrix of parameter estimates (row 1) and corresponding standard errors (row 2). Parameters estimated are 1. \eqn{\lambda}, the incidence rate; 
2. \eqn{\phi}, the proportion of the population uninfected at the time of cross-sectional sample (1-prevalence); 3. \eqn{\mu}, the mean time spent in the recent infection state; 4. 
\eqn{p}, proportion of infected subjects who will eventually test positive on the less-sensitive test. }
\item{incidence.ci }{A \code{100*(1-alpha)\%} confidence interval for \eqn{\lambda}, calculated by transforming and back-transforming the \code{100*(1-alpha)\%} confidence interval for 
\eqn{\log(\lambda)} via the delta method, to ensure that the interval contains only positive values. This interval is symmetric about  
\eqn{\log(\hat{\lambda})}, and not about \eqn{\hat{\lambda}}.}
}
\references{
  Wang R and Lagakos S. (2009).
  On the Use of Adjusted Cross-Sectional Estimators of HIV Incidence.
  \emph{J Acquir Immune Defic Syndr}, 52(5):538-547.

  Wang R and Lagakos S. (2010).
  Augmented Cross-Sectional Prevalence Testing for Estimating HIV Incidence.
  \emph{Biometrics}, 66:864-874.

  Claggett B, Lagakos S, Wang R. (2011).
  Augmented Cross-Sectional Studies with Abbreviated Follow-up for Estimating HIV Incidence.
  \emph{Biometrics}, DOI:10.1111/j.1541-0420.2011.01632.x.
}
\author{
Brian Claggett <bclagget@hsph.harvard.edu>,\cr
Weiliang Qiu <stwxq@channing.harvard.edu>,\cr
Rui Wang <rwang@hsph.harvard.edu>
}
\examples{
####Example with Simulated data

set.seed(1234)
N = 6000
phi = 0.85
lambda = 0.02
Mu = 0.5
p = 0.95
W = 52  ###max followup time
interval = 2 ### interval between follow-up visits, in weeks

pt = (p * phi * lambda * Mu) / (p * phi * lambda * Mu + (1 - p) * (1 - phi))


Ns =  rmultinom(1, N, 
          c(phi, phi * p * lambda * Mu, (1 - p) * (1 - phi), 
            p * (1 - phi - phi * lambda * Mu)))

t1 = runif(500000,  -5, 0)
t2 = t1 + rweibull(500000, 2.7, Mu / gamma(1 + 1 / 2.7))
frd = t2[t2 > 0]
X = sample(frd, Ns[2])

### interval censoring ####
a = floor(X * (52 / interval)) * interval
b = ceiling(X * (52 / interval)) * interval
ints = cbind(a, b)
sortints = ints[order(ints[, 1]), ]

N2 = Ns[2] + Ns[3]
n1t = sum(sortints[, 2] <= W)

################################################
## p=1; in the absence of long-term non-progressors ##
csiaug.est(N1 = Ns[1], n1 = N2, n0 = 0, N3 = Ns[4], 
    ftime = rbind( sortints[1:n1t, ], 
                   t(matrix(rep(c(W, 999), N2 - n1t), nrow = 2)) ) / 52,
    alpha = 0.05)

################################################
## p<1; in the presence of long-term non-progressors ##
csiaug.est(N1 = Ns[1], n1 = n1t, n0 = N2 - n1t, N3 = Ns[4], 
    ftime = sortints[1:n1t, ] / 52,
    alpha = 0.05)

}
% Add one or more standard keywords,  see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }


