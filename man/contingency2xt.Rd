\name{contingency2xt}
\alias{contingency2xt}
\title{
Kruskal-Wallis Test for the 2 x t Contingency Table 
}
\description{
This function uses the Kruskal-Wallis criterion to test
the hypothesis of no association between the counts 
for two responses
"A" and "B" across t categories.
}
\usage{
contingency2xt(Avec, Bvec, 
	method = c("asymptotic", "simulated", "exact"), 
	dist = FALSE, tab0 = TRUE, Nsim = 1e+06)
}
\arguments{
  	\item{Avec}{
		vector of length \eqn{t} giving the counts \eqn{A_1,\ldots, A_t}
		for response "A" according to \eqn{t} categories. 
		\eqn{m = A_1 + \ldots + A_t}.
	}
  	\item{Bvec}{
		vector of length \eqn{t} giving the counts \eqn{B_1,\ldots, B_t}
		for response "B" according to \eqn{t} categories. 
		\eqn{n = B_1 + \ldots + B_t = N-m}.
	}
    \item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic chi-square approximation 
		with \eqn{t-1} degrees of freedom to approximate the \eqn{P}-value. 
		This calculation is always done.

		\code{"simulated"} uses \code{Nsim} simulated counts for \code{Avec} and 
		\code{Bvec} with the observed marginal totals, \code{m, n, d = Avec+Bvec}, 
		to estimate the \eqn{P}-value.

		\code{"exact"} enumerates all counts for \code{Avec} and \code{Bvec} with 
		the observed marginal totals to get an exact \eqn{P}-value. It is used only 
		when \code{Nsim} is at least as large as the number \code{choose(m+t-1,t-1)}
		of full enumerations. 
		Otherwise, \code{method} reverts to \code{"simulated"} using the given \code{Nsim}.
    }
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{dist = TRUE}, the distribution of the
		simulated or fully enumerated Kruskal-Wallis test statistics is
		returned as \code{null.dist}, if \code{dist = FALSE} the value
		of \code{null.dist} is \code{NULL}.
		The coice \code{dist = TRUE} also limits \code{Nsim <- min(Nsim,1e8)}.
	}
	\item{tab0}{\code{TRUE} (default) or \code{FALSE}. If \code{tab0 = TRUE}, the null distribution
		is returned in 2 column matrix form when 
		\code{method = "simulated"}. When \code{tab0 = FALSE} the simulated null distribution 
		is returned as a vector of all simulated values of the test statistic.
	}
	\item{Nsim}{\code{=10000} (default), number of simulated \code{Avec} splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method =}
		\code{"simulated"}, as previously explained.
	}

}
\details{
For this data scenario the Kruskal-Wallis criterion is
\deqn{K.star = \frac{N(N-1)}{mn}(\sum\frac{A_i^2}{d_i}-\frac{m^2}{N})}{K.star = 
N(N-1)/(mn) (\sum A_i^2/d_i-m^2/N)}
with \eqn{d_i=A_i+B_i}, treating "A" responses
as 1 and "B" responses as 2, and using midranks as explained in Lehmann (2006), Chapter 5.3.

For small sample sizes exact null distribution
calculations are possible, based on Algorithm C (Chase's sequence) in Knuth (2011),
which allows the enumeration of all possible splits of \eqn{m} into counts
\eqn{A_1,\ldots, A_t} such that 
\eqn{m = A_1 + \ldots + A_t},
followed by the calculation of the statistic 
\eqn{K.star} for each such split. 
Simulation of \eqn{A_1,\ldots, A_t} uses the probability model (5.35) in Lehmann (2006)
to successively generate hypergeometric counts \eqn{A_1,\ldots, A_t}.
Both these processes, enumeration and simulation, are done in C.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"2 x t Contingency Table"}}
\item{t}{number of classification categories}
\item{KW.cont}{2 (3) vector giving the observed KW statistic, its asymptotic 
\eqn{P}-value (and simulated or exact \eqn{P}-value) }
\item{null.dist}{simulated or enumerated null distribution 
of the test statistic. It is given as an \code{M} by 2 matrix,
where the first column (named \code{KW}) gives the \code{M} unique ordered 
values of the Kruskal-Wallis 
statistic and the second column (named \code{prob}) gives the corresponding (simulated or exact)
probabilities. 

This format of \code{null.dist} is returned when \code{method = "exact"} 
and \code{dist} \code{= TRUE} or when \code{method =}\code{ "simulated"} 
and \code{dist = TRUE} and \code{tab0} \code{= TRUE} are specified.

For \code{method =} \code{"simulated"}, \code{dist = TRUE}, and 
\code{tab0 = FALSE} the null distribution \code{null.dist} is returned as the vector of
all simulated test statistic values. This is used in \code{\link{contingency2xt.comb}}
in the simulation mode.

\code{null.dist = NULL} is returned 
when \code{dist = FALSE} or when \code{method =}
\code{"asymptotic"}. }
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations.}
}

\section{warning}{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. In most cases
\code{dist = TRUE} should not be used, i.e.,  
when the returned distribution objects 
become too large for R's work space.}

\references{
Knuth, D.E. (2011), \emph{The Art of Computer Programming, Volume 4A 
Combinatorial Algorithms Part 1}, Addison-Wesley


Kruskal, W.H. (1952), A Nonparametric Test for the Several Sample Problem,
\emph{The Annals of Mathematical Statistics},
\bold{Vol 23, No. 4}, 525-540

Kruskal, W.H. and Wallis, W.A. (1952), Use of Ranks in One-Criterion Variance Analysis,
\emph{Journal of the American Statistical Association}, 
\bold{Vol 47, No. 260}, 583--621. 

Lehmann, E.L. (2006), \emph{Nonparametrics, Statistical Methods Based on Ranks},
Revised First Edition, 
Springer, New York.

}




\examples{
contingency2xt(c(25,15,20),c(16,6,18),method="exact",dist=FALSE,
	tab0=TRUE,Nsim=1e3)
}

\keyword{nonparametric}
\keyword{htest}

