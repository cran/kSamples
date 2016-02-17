\name{pp.kSamples}
\alias{pp.kSamples}
\title{
Upper Tail Probability Plots for Objects of Class kSamples
}
\description{
This function plots upper tail probabilities of the limiting distribution 
against the corresponding exact or simulated probabilities, both on a log-scale.
}
\usage{
pp.kSamples(x) 
}
\arguments{
  	\item{x}{an object of class \code{kSamples}}
	}

\details{Objects of class \code{kSamples} are produced by any of the following functions

\code{\link{ad.test}} Anderson-Darling k-sample test.

\code{\link{ad.test.combined}} Combined Anderson-Darling k-sample tests.

\code{\link{qn.test}} \eqn{QN} rank scores test.

\code{\link{qn.test.combined}} Combined \eqn{QN} rank scores tests.

\code{\link{contingency2xt}} test for \eqn{2 * t} contingency table.

\code{\link{contingency2xt.comb}} test for the combination of \eqn{2 * t} contingency tables.

\code{\link{jt.test}} Jonckheere-Terpstra test.

\code{\link{Steel.test}} Steel test. This will work only for alternative = "greater" or "two-sided".
The approximation quality for "less" is the same as for "greater".

The command \code{pp.kSamples(x)} for an object of class \code{kSamples}
will only produce a plot when the object \code{x} contains 
non-NULL entries for the null distribution. The purpose of this function is to give the user
a sense of the asymptotic distribution accuracy.
}


\seealso{
\code{\link{ad.test}},
\code{\link{ad.test.combined}},
\code{\link{qn.test}},
\code{\link{qn.test.combined}},

\code{\link{contingency2xt}},
\code{\link{contingency2xt.comb}}
\code{\link{jt.test}}
\code{\link{Steel.test}}
}

\examples{
qn.out <- qn.test(c(1,3,7,2,9),c(1,4,6,11,2),test="KW",
		method="simulated",dist=TRUE,Nsim=1000) 
pp.kSamples(qn.out)
}

\keyword{nonparametric}
\keyword{htest}

