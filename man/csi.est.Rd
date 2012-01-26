\name{csi.est}
\alias{csi.est}
\title{
Cross-sectional incidence estimation for HIV surveillance.
}
\description{
Estimating incidence rate from assay-based cross-sectional studies, allowing for a subpopulation that remains negative on the less-sensitive assay indefinitely.}
\usage{
csi.est(N1, N2, N3, mu = 0.5, p = 1, se.mu = 0, se.p = 0, alpha = 0.05)
}
\arguments{
  \item{N1}{
Number of subjects testing negative on the sensitive test (e.g. ELISA).
}
  \item{N2}{
Number of subjects testing positive on the sensitive test (e.g. ELISA) but
Testing negative on the less-sensitive test (e.g. BED).
}
  \item{N3}{
Number of subjects testing positive on the less-sensitive test (e.g. BED).
%%     ~~Describe \code{n3} here~~
}
  \item{mu}{
Mean window period for the subpopulation of infected individuals
that would eventually test positive on the less-sensitive test, expressed on the same time scale for which incidence rate is desired. 
}
  \item{p}{
Proportion of subjects who would become testing positive at some point
after seroconversion. \code{p=1} corresponds to the setting where all infected subjects will test positive on the less-sensitive test. 
}
  \item{se.mu}{
standard error of an externally-estimated \eqn{\hat{\mu}}. 
}
  \item{se.p}{
standard error of an externally-estimated \eqn{\hat{\mu}}. 
}
  \item{alpha}{
specifys the coverage level (\code{100*(1-alpha)\%}) of resulting confidence interval for incidence rate.
}
}
\details{
Confidence interval is obtained via the delta method, and is symmetric on the log scale.
}
\value{
A list with 3 elements:
\item{inc.est }{estimated incidence rate}
\item{se.est }{standard error of the estimated incidence rate}
\item{incidence.ci}{\code{100*(1-alpha)\%}  confidence interval of the incidence rate}
}
\references{
  Wang R and Lagakos S. (2009).
  On the Use of Adjusted Cross-Sectional Estimators of HIV Incidence.
  \emph{J Acquir Immune Defic Syndr}, 52(5):538-547.
}
\author{
Brian Claggett <bclagget@hsph.harvard.edu>,\cr
Weiliang Qiu <stwxq@channing.harvard.edu>,\cr
Rui Wang <rwang@hsph.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
csi.est(N1 = 5126, N2 = 77, N3 = 797, mu = 0.5, p = 1, 
    se.mu = 0, se.p = 0, alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

