\name{qn.test.combined}
\alias{qn.test.combined}
\title{
Combined Rank Score k-Sample Tests
}
\description{
This function combines several independent rank score \eqn{k}-sample tests
into one overall test of the hypothesis that the independent samples within 
each block come from a common unspecified distribution, while the common
distributions may vary from block to block.  
}
\usage{
qn.test.combined(\dots, data = NULL, test = c("KW", "vdW", "NS"),
	method = c("asymptotic", "simulated", "exact"),
	dist = FALSE, Nsim = 10000)
}
\arguments{
  	\item{\dots}{
		Either a sequence of several lists, say \eqn{L_1, \ldots, L_M} (\eqn{M > 1})
		where list \eqn{L_i} contains \eqn{k_i > 1}  sample vectors of respective 
		sizes \eqn{n_{i1}, \ldots, n_{ik_i}}, 
		where \eqn{n_{ij} > 4} is recommended
		for reasonable asymptotic \eqn{P}-value calculation. 
		\eqn{N_i=n_{i1}+\ldots+n_{ik_i}}
		is the pooled sample size for block \eqn{i},

		or a list of such lists,

		or a formula, like y ~ g | b, where y is a numeric response vector,
		g is a factor with levels indicating different treatments and
		b is a factor indicating different blocks; y, g, b have same length. 
		y is split separately for each block level into separate samples 
		according to the g levels. The same g level may occur in different 
		blocks. The variable names may correspond to variables in an optionally 
		supplied data frame via the data = argument.
	}
    	\item{data}{= an optional data frame providing the variables in formula input
		}
	\item{test}{= \code{c("KW", "vdW", "NS")}, 

		where \code{"KW"} uses scores \code{1:N} (Kruskal-Wallis test)

		\code{"vdW"} uses van der Waerden scores, \code{qnorm( (1:N) / (N+1) )}

		\code{"NS"} uses normal scores, i.e., expected values of 
		standard normal order statistics,
		invoking function \code{normOrder} of \code{package SuppDists (>=1.1-9.4)}

		For the above scores \eqn{N} changes from block to block and represents the total 
		pooled sample size \eqn{N_i}{N.i} for block \eqn{i}. 
	}
    \item{method}{=\code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic chi-square approximation for the  \eqn{P}-value.
		The adequacy of asymptotic \eqn{P}-values for use with moderate sample sizes
		may be checked with \code{method = "simulated"}.

		\code{"simulated"} uses simulation to get \code{Nsim} simulated \eqn{QN} statistics for each
		block of samples, adding them component wise across blocks to get \code{Nsim} 
		combined values, and compares these with the observed combined value to 
		get the estimated \eqn{P}-value.

		\code{"exact"} uses full enumeration of the test statistic value 
		for all sample splits of the pooled samples within each block.
		The test statistic vectors for each block are added 
		(each component against each component, as in the R \code{outer(x,y,} \code{"+")} command) 
		to get the convolution enumeration for the combined test statistic.
		This "addition" is done one block at a time. 
		It is possible only for small problems, and is attempted only when \code{Nsim}
		is at least the (conservatively maximal) length 
		\deqn{\frac{N_1!}{n_{11}!\ldots n_{1k_1}!}\times\ldots\times \frac{N_M!}{n_{M1}!\ldots n_{Mk_M}!}}
		of the final distribution vector, were \eqn{N_i = n_{i1}+\ldots+n_{ik_i}} 
		is the sample size of the pooled samples for the i-th block. Otherwise, it reverts to the 
		simulation method using the provided \code{Nsim}.
	}
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated convolution vector \code{null.dist} is returned for the 
		\eqn{QN} test statistic. 

		Otherwise, \code{NULL} is returned.
	}
	\item{Nsim}{\code{= 10000} (default), number of simulation splits to use within 
		each block of samples. It is only used when \code{method =} \code{"simulated"}
		or when \code{method =} \code{"exact"} reverts to \code{method =} \code{"simulated"}, 
		as previously explained.
		Simulations are independent across blocks, using \code{Nsim} for each block.
	}

}
\details{
The rank score \eqn{QN} criterion \eqn{QN_i}{QN.i} for the \eqn{i}-th block of 
\eqn{k_i}{k.i} samples, 
is used to test the hypothesis that the samples in the \eqn{i}-th block all come 
from the same but unspecified continuous distribution function \eqn{F_i(x)}{F.i(x)}.
See \code{\link{qn.test}} for the definition of the \eqn{QN} criterion
and the exact calculation of its null distribution.

The combined \eqn{QN} criterion \eqn{QN_{\rm comb} = QN_1 + \ldots + QN_M}{%
QN.comb = QN.1 + \ldots + QN.M}
is used to simultaneously test whether the samples 
in block i come from the same continuous distribution function \eqn{F_i(x)}{F.i(x)}. 
However, the unspecified common distribution function \eqn{F_i(x)}{F.i(x)} may change 
from block to block. 

The \eqn{k} for each block of \eqn{k}
independent samples may change from block to block.

The asymptotic approximating chi-square distribution has 
\eqn{f = (k_1-1)+\ldots+(k_M-1)}{f = (k.1-1)+\ldots+(k.M-1)} degrees of freedom.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

The continuity assumption can be dispensed with if we deal with 
independent random samples, or if randomization was used in allocating
subjects to samples or treatments, independently from block to block, and if
the asymptotic, simulated or exact \eqn{P}-values are viewed conditionally, given the tie patterns
within each block. Under such randomization any conclusions 
are valid only with respect to the blocks of subjects that were randomly allocated.
In case of ties the average rank scores are used across tied observations within each block.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Kruskal-Wallis"}, \code{"van der Waerden scores"}, or

\code{"normal scores"}}
\item{M}{number of blocks of samples being compared}
\item{n.samples}{list of \code{M} vectors, each vector giving the sample sizes for 
each block of samples being compared}
\item{nt}{vector of length \code{M} of total sample sizes involved in each of the 
\code{M} comparisons of \eqn{k_i}{k.i} samples, respectively}
\item{n.ties}{vector giving the number of ties in each the \code{M}
comparison blocks}
\item{qn.list}{list of \code{M} matrices giving the \code{qn} results 
from \code{qn.test}, applied to the samples in each of
the \code{M} blocks}
\item{qn.c}{2 (or 3) vector containing the observed 
\eqn{QN_{\rm comb}}{QN.comb}, asymptotic \eqn{P}-value, 
(simulated or exact \eqn{P}-value).}
\item{warning}{logical indicator, \code{warning = TRUE} when at least one 
\eqn{n_{ij} < 5}{n.ij < 5}.}
\item{null.dist}{simulated or enumerated null distribution of the 
\eqn{QN_{\rm comb}}{QN.comb} statistic.
It is \code{NULL} when \code{dist = FALSE} or when \code{method = "asymptotic"}.}
\item{method}{The \code{method} used.}
\item{Nsim}{The number of simulations used for each block of samples.}
}
\references{
Lehmann, E.L. (2006), \emph{Nonparametric, Statistical Methods Based on Ranks}, Springer Verlag, New York. Ch. 6, Sec. 5D.
}

\note{
These tests are useful in analyzing treatment effects of shift nature in randomized 
(incomplete) block experiments.
}
\seealso{
\code{\link{qn.test}}
}
\examples{
## Create two lists of sample vectors.
x1 <- list( c(1, 3, 2, 5, 7), c(2, 8, 1, 6, 9, 4), c(12, 5, 7, 9, 11) )
x2 <- list( c(51, 43, 31, 53, 21, 75), c(23, 45, 61, 17, 60) )
# and a corresponding data frame datx1x2
x1x2 <- c(unlist(x1),unlist(x2))
gx1x2 <- as.factor(c(rep(1,5),rep(2,6),rep(3,5),rep(1,6),rep(2,5)))
bx1x2 <- as.factor(c(rep(1,16),rep(2,11)))
datx1x2 <- data.frame(A = x1x2, G = gx1x2, B = bx1x2)

## Run qn.test.combined.
set.seed(2627)
qn.test.combined(x1, x2, method = "simulated", Nsim = 1000) 
# or with same seed
# qn.test.combined(list(x1, x2), method = "simulated", Nsim = 1000)
# or qn.test.combined(A~G|B,data=datx1x2,method="simulated",Nsim=1000)
}



\keyword{nonparametric}
\keyword{htest}
\keyword{design}
