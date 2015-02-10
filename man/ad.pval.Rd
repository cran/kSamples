\name{ad.pval}
\alias{ad.pval}
\title{
\eqn{P}-Value for the Asymptotic Anderson-Darling Test Distribution
}
\description{
This function computes upper tail probabilities for the limiting distribution of the 
standardized Anderson-Darling test statistic.
}
\usage{
ad.pval(tx,m,version=1) 
}
\arguments{
  	\item{tx}{
		a vector of desired thresholds \eqn{\ge 0}
	}
    \item{m}{ 
		The degrees of freedom for the asymptotic standardized Anderson-Darling test statistic
    }
	\item{version}{\code{= 1} (default) if \eqn{P}-value for version 1 of the 
		test statistic is desired, otherwise the version 2 \eqn{P}-value is calculated.
	}
}
\details{
Extensive simulations (sampling from a common continuous distribution)
were used to extend the range of the asymptotic 
\eqn{P}-value calculation from the original \eqn{[.01,.25]} in Table 1 of the reference paper
to 36 quantiles corresponding to \eqn{P}  = .00001, .00005, .0001, .0005, .001, .005, .01, .025, 
.05, .075, .1, .2, .3, .4, .5, .6, .7, .8, .9, .925, .95, .975, .99, .9925, .995, .9975, .999,
.99925, .9995, .99975, .9999, .999925, .99995, .999975, .99999. Note that the entries of the original Table 1
were obtained by using the first 4 moments of the asymptotic distribution and a
Pearson curve approximation.

Using \code{ad.test},
1 million replications of the standardized \eqn{AD} statistics with sample sizes
\eqn{n_i=500}{n.i=500}, \eqn{i=1,\ldots,k} were run for \eqn{k=2,3,4,5,7} (\eqn{k=2} was done twice).
These values of \eqn{k} correspond to degrees of freedom
\eqn{m=k-1=1,2,3,4,6} in the asymptotic distribution. The random variable described by this 
distribution is denoted by \eqn{T_m}{T.m}.
The actual variances (for \eqn{n_i=500}{n.i=500}) agreed fairly well with the asymptotic variances.

Using the convolution nature of the asymptotic distribution, the performed simulations
were exploited to result in an effective simulation of 2 million cases, except for 
\eqn{k=11}, i.e., \eqn{m=k-1=10}, for which the asymptotic distribution of 
\eqn{T_{10}}{T.10}
was approximated by the sum of the \eqn{AD} statistics for \eqn{k=7} and \eqn{k=5}, 
for just the 1 million cases run for each \eqn{k}. 

The interpolation of tail 
probabilities \eqn{P}
for any desired \eqn{k} is done in two stages. First, a spline in \eqn{1/\sqrt m} is
fitted to each of the 36 quantiles obtained for \eqn{m=1,2,3,4,6,8,10,\infty} to obtain 
the corresponding interpolated quantiles for the \eqn{m} in question. 

Then a spline is fitted
to the \eqn{\log((1-P)/P)}{log((1-P)/P)} as a function of these 36 interpolated quantiles. This latter
spline is used to determine the tail probabilities \eqn{P} for the 
specified threshold \code{tx}, corresponding to either \eqn{AD} 
statistic version. The above procedure is based on simulations for either version
of the test statistic,
appealing to the same limiting distribution.
}
\value{
a vector of upper tail probabilities corresponding to \code{tx}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, 
\emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}



\seealso{
\code{\link{ad.test}}
\code{\link{ad.test.combined}}
}

\examples{
ad.pval(tx=c(3.124,5.65),m=2,version=1)
ad.pval(tx=c(3.124,5.65),m=2,version=2)
}

\keyword{nonparametric}
\keyword{htest}

