\name{contingency2xt.comb}
\alias{contingency2xt.comb}
\title{
Combined Kruskal-Wallis Tests for the 2 x t Contingency Tables
}
\description{
This function uses the Kruskal-Wallis criterion to test
the hypothesis of no association between the counts 
for two responses
"A" and "B" across t categories and across \eqn{M} blocks.
}
\usage{
contingency2xt.comb(\dots, 
method = c("asymptotic", "simulated", "exact"), 
dist = FALSE, Nsim = 10000)
}
\arguments{
  	\item{\dots}{ 
		Either several lists \eqn{L_1,\ldots,L_M}{L.1,\ldots,L.M}, each
		of two equal length vectors \eqn{A_i}{A.i} and
		\eqn{B_i}{B.i}, \eqn{i=1,\ldots,M}, of counts \eqn{\ge 0}, 
		where the common length \eqn{t_i}{t.i} of \eqn{A_i}{A.i} and
		\eqn{B_i}{B.i} may vary from list to list

		or a list of \code{M} such lists
	}
    \item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic chi-square approximation 
		with \eqn{(t_1-1)+\ldots+(t_M-1)}{(t.1-1)+\ldots+(t.M-1)} degrees of freedom
		to approximate the \eqn{P}-value, This calculation is always done.

		\code{"simulated"} uses \code{Nsim} simulated counts for the two vectors
		\eqn{A_i}{A.i} and \eqn{B_i}{B.i} in list \eqn{L_i}{L.i}, 
		with the observed marginal totals, \eqn{m_i=\sum A_i}{m.i=\sum A.i}, 
		\eqn{n_i = \sum B_i}{n.i=\sum B.i}, \eqn{d_i = A_i+B_i}{d.i=A.i+B.i}. 
		It does this independently from list to list using the same \code{Nsim} each time, 
		adding the resulting Kruskal-Wallis criteria across lists
		to get \code{Nsim} such summed values to estimate the \eqn{P}-value.

		\code{"exact"} enumerates all counts for \eqn{A_i}{A.i} and \eqn{B_i}{B.i} with 
		the respective observed marginal totals to get an exact distribution for each list. 
		These distributions are then convolved to obtain the \eqn{P}-value.
		It is used only when \code{Nsim} is at least as large as the product across blocks
		of the number \code{choose(m+t-1,t-1)} of full enumerations per block.
		Otherwise, \code{method} reverts to \code{"simulated"} using the given \code{Nsim}.
    }
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated null distribution \code{null.dist} is returned
		for the Kruskal-Wallis test statistic. Otherwise \code{null.dist = NULL} is returned.
	}
	\item{Nsim}{\code{=10000} (default), number of simulated \eqn{A_i}{A.i} splits to use per block.
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method = "simulated"}, as previously explained.
	}

}
\details{
For details on the calculation of the Kruskal-Wallis criterion and its exact or simulated
distribution for each block see \code{\link{contingency2xt}}.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Combined 2 x t Contingency Tables"}}
\item{t}{vector giving the number of classification categories per block}
\item{M}{number of blocked tables}
\item{kw.list}{a list of the \code{KW.cont} output componenents from 
\code{\link{contingency2xt}} for each of the blocks}
\item{null.dist}{simulated or enumerated null distribution 
of the combined test statistic. It is given as an \code{L} by 2 matrix,
where the first column (named \code{KW}) gives the \code{L} unique ordered 
values of the combined Kruskal-Wallis 
statistic and the second column (named \code{prob}) gives the corresponding (simulated or exact)
probabilities.

\code{null.dist = NULL} is returned when \code{dist = FALSE} or when 
\code{method =} \code{"asymptotic"}.}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations.}
}
\note{
The required level for \code{Nsim} in order for \code{method = "exact"}
to be carried out, is conservative, but there is no transparent way to get a 
better estimate. The actual dimension \code{L} of the realized \code{null.dist}
will typically be much smaller, since the distribution is compacted to
its unique support values.
}

\section{warning}{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. In most cases
\code{dist = TRUE} should not be used, i.e.,  
when the returned distribution objects 
become too large for R's work space.}



\examples{
out <- contingency2xt.comb(list(c(25,15,20),c(16,6,18)),
list(c(12,4,5),c(13,8,9)),method = "simulated", dist=FALSE, Nsim=1e3)
}

\keyword{nonparametric}
\keyword{htest}

