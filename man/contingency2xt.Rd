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
contingency2xt(Avec, Bvec, method = 
	c("asymptotic", "simulated", "exact"), dist = FALSE, 
	tab0 = TRUE, Nsim = 1e+06)
}
\arguments{
  	\item{Avec}{
		vector of length \eqn{t} giving the counts \eqn{A_1,\ldots, A_t}{A.1,\ldots,A.t} 
		for response "A" according to \eqn{t} categories. 
		\eqn{m = A_1 + \ldots + A_t}{m = A.1 + \ldots + A.t}.
	}
  	\item{Bvec}{
		vector of length \eqn{t} giving the counts \eqn{B_1,\ldots, B_t}{B.1,\ldots,B.t}
		for response "B" according to \eqn{t} categories. 
		\eqn{n = B_1 + \ldots + B_t = N-m}{n = B.1 + \ldots + B.t}.
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
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, the distribution of the
		simulated or fully enumerated Kruskal-Wallis test statistics is
		returned as \code{null.dist}.
	}
	\item{tab0}{\code{TRUE} (default) or \code{FALSE}. If \code{TRUE}, the null distribution
		is returned in 2 column matrix form when 
		\code{method = "simulated"}. When \code{tab0 = FALSE} the simulated null distribution 
		is returned as a vector of all simulated values of the test statistic.
	}
	\item{Nsim}{\code{=10000} (default), number of simulated \code{Avec} splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method = "simulated"}, as previously explained.
	}

}
\details{
For this data scenario the Kruskal-Wallis criterion is
\deqn{\bar{K}^\ast = \frac{N(N-1)}{mn}\left(\sum\frac{A_i^2}{d_i}-\frac{m^2}{N}\right)}{%
KW = [ N(N-1)/(mn) ] * [ \sum (A.i^2)/d.i - (m^2)/N ]}
with \eqn{d_i=A_i+B_i}{d.i=A.i+B.i}, treating "A" responses
as 1 and "B" responses as 2, and using midranks as explained in Lehmann (2006), Chapter 5.3.

For small sample sizes exact null distribution
calculations are possible, based on Algorithm C (Chase's sequence) in Knuth (2011),
which allows the enumeration of all possible splits of \eqn{m} into counts
\eqn{A_1,\ldots, A_t}{A.1,\ldots,A.t} such that 
\eqn{m = A_1 + \ldots + A_t}{m = A.1 + \ldots + A.t},
followed by the calculation of the statistic 
\eqn{\bar{K}^\ast}{KW} for each such split. 
Simulation of \eqn{A_1,\ldots, A_t}{A.1,\ldots, A.t} uses the probability model (5.35) in Lehmann (2006)
\deqn{P(A_1=a_1,\ldots,A_t=a_t)=\frac{{d_1\choose a_1}\cdots {d_t\choose a_t}}{{N\choose m}}}{%
P(A.1=a.1,\ldots,A.t=a.t)=choose(d.1,a.1) * \ldots * choose(d.t,a.t) / choose(N,m)}
\deqn{
=\frac{{d_1\choose a_1}{N-d_1\choose m-a_1}}{{N\choose m}}
\frac{{d_2\choose a_2}{N-d_1-d_2\choose m-a_1-a_2}}{{N-d_1\choose m-a_1}}
\frac{{d_3\choose a_3}{N-d_1-d_2-d_3\choose m-a_1-a_2-a_3}}{{N-d_1-d_2\choose
 m-a_1-a_2}}
 \ldots
\frac{{d_t\choose a_t}{N-d_1-d_2-d_3-\dots-d_t\choose m-a_1-a_2-a_3-
\ldots-a_t}}{{N-d_1-d_2-\ldots-d_{t-1}\choose
 m-a_1-a_2-\ldots-a_{t-1}}}}{%
= [ choose(d.1,a.1)*choose(N-d.1,m-a.1) / choose(N,m) ]}
\deqn{}{ *
[ choose(d.2,a.2)*choose(N-d.1-d.2,m-a.1-a.2) / choose(N-d.1, m-a.1) ] }
\deqn{}{ *
[ choose(d.3,a.3)*choose(N-d.1-d.2-d.3,m-a.1-a.2-a.3) / choose(N-d.1-d.2, m-a.1-a.2) ]
* \ldots 
}
to successively generate hypergeometric counts \eqn{A_1,\ldots, A_t}{A.1,\ldots,A.t}.
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

