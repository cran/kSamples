\name{ad.test.combined}
\alias{ad.test.combined}
\title{
Combined Anderson-Darling k-Sample Tests
}
\description{
This function combines several independent Anderson-Darling \eqn{k}-sample tests
into one overall test of the hypothesis that the independent samples within 
each block come from a common unspecified distribution, while the common
distributions may vary from block to block.  Both versions of the 
Anderson-Darling test statistic are provided.
}
\usage{
ad.test.combined(\dots, data = NULL,
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
		b is a factor indicating different blocks; y, g, b are or equal length. 
		y is split separately for each block level into separate samples 
		according to the g levels. The same g level may occur in different blocks. 
		The variable names may correspond to variables in an optionally supplied 
		data frame via the data = argument,
	}
    \item{data}{= an optional data frame providing the variables in formula input
		}
    \item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic \eqn{P}-value approximation, reasonable 
		for P in [0.00001, .99999], linearly extrapolated via 
		\eqn{\log(P/(1-P))}{log(P/(1-P))} outside 
		that range. See \code{\link{ad.pval}} for details. 		
		The adequacy of the asymptotic \eqn{P}-value calculation may be checked using 
		\code{\link{pp.kSamples}}.


		\code{"simulated"} uses simulation to get \code{Nsim} simulated \eqn{AD} statistics 
		for each block of samples, adding them across blocks component wise to get \code{Nsim} 
		combined values. These are compared with the observed combined value to obtain the 
		estimated	\eqn{P}-value.

		\code{"exact"} uses full enumeration of the test statistic values 
		for all sample splits of the pooled samples within each block.
		The test statistic vectors for the first 2 blocks are added 
		(each component against each component, as in the R \code{outer(x,y,} \code{"+")} command) 
		to get the convolution enumeration for the combined test statistic. The resulting
		vector is convoluted against the next block vector in the same fashion, and so on.
		It is possible only for small problems, and is attempted only when \code{Nsim}
		is at least the (conservatively maximal) length 
		\deqn{\frac{N_1!}{n_{11}!\ldots n_{1k_1}!}\times\ldots\times 
			\frac{N_M!}{n_{M1}!\ldots n_{Mk_M}!}}
		of the final distribution vector. Otherwise, it reverts to the 
		simulation method using the provided \code{Nsim}.
    	}
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated convolution vectors 
		\code{null.dist1} and \code{null.dist2} are returned for the respective
		test statistic versions. Otherwise, \code{NULL} is returned for each.
	}
	\item{Nsim}{\code{= 10000} (default), number of simulation splits to use within 
		each block of samples. It is only used when \code{method = "simulated"}
		or when \code{method =} \code{"exact"} reverts to \code{method = "simulated"}, 
		as previously explained. Simulations are independent across blocks, 
		using \code{Nsim} for each block. \code{Nsim} is limited by \code{1e7}.
	}

}
\details{
If \eqn{AD_i} is the Anderson-Darling criterion for the i-th block of 
\eqn{k_i} samples, 
its standardized test statistic is 
\eqn{T_i = (AD_i - \mu_i)/\sigma_i}, with 
\eqn{\mu_i} and
\eqn{\sigma_i} representing mean and standard deviation of 
\eqn{AD_i}. This statistic 
is used to test the hypothesis that the samples in the i-th block all come 
from the same but unspecified continuous distribution function \eqn{F_i(x)}.

The combined Anderson-Darling criterion is 
\eqn{AD_{comb}=AD_1 + \ldots + AD_M} and 
\eqn{T_{comb} = } \eqn{(AD_{comb} - \mu_c)/\sigma_c} is the standardized form, 
where \eqn{\mu_c=\mu_1+\ldots+\mu_M} and \eqn{\sigma_c =
\sqrt{\sigma_1^2 +\ldots+\sigma_M^2}} 
represent the mean and standard deviation of \eqn{AD_{comb}}.
The statistic \eqn{T_{comb}} is used to simultaneously 
test whether the samples 
in each block come from the same continuous distribution function 
\eqn{F_i(x), i=1,\ldots,M}. 
The unspecified common distribution function \eqn{F_i(x)} may change 
from block to block. According to the reference article, two versions
of the test statistic and its corresponding combinations are provided.

The \eqn{k_i} for each block of \eqn{k_i}
independent samples may change from block to block.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

The continuity assumption can be dispensed with if we deal with 
independent random samples, or if randomization was used in allocating
subjects to samples or treatments, independently from block to block, and if we view
the simulated or exact \eqn{P}-values conditionally, given the tie patterns
within each block. Of course, under such randomization any conclusions 
are valid only with respect to the blocks of subjects that were randomly allocated.
The asymptotic \eqn{P}-value calculation assumes distribution continuity. No adjustment 
for lack thereof is known at this point. The same comment holds for the means
and standard deviations of respective statistics.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\eqn{=} \code{"Anderson-Darling"}}
\item{M}{number of blocks of samples being compared}
\item{n.samples}{list of \code{M} vectors, each vector giving the sample sizes for 
each block of samples being compared}
\item{nt}{\eqn{= (N_1,\ldots,N_M)}}
\item{n.ties}{vector giving the number of ties in each the \code{M}
comparison blocks}
\item{ad.list}{list of \code{M} matrices giving the \code{ad} results 
for \code{ad.test} applied to the samples in each of
the \code{M} blocks}
\item{mu}{vector of means of the \eqn{AD} statistic for the  \code{M} blocks}
\item{sig}{vector of standard deviations of the \eqn{AD} statistic for the  \code{M} blocks}
\item{ad.c}{2 x 3 (2 x 4) matrix containing 
\eqn{AD_{comb}, T_{comb}}, asymptotic \eqn{P}-value, 
(simulated or exact \eqn{P}-value), for each version of the combined test statistic,
version 1 in row 1 and version 2 in row 2}
\item{mu.c}{mean of \eqn{AD_{comb}}}
\item{sig.c}{standard deviation of \eqn{AD_{comb}}}
\item{warning}{logical indicator, \code{warning = TRUE} when at least one 
\eqn{n_{ij} < 5}}
\item{null.dist1}{simulated or enumerated null distribution of version 1 
of \eqn{AD_{comb}}}
\item{null.dist2}{simulated or enumerated null distribution of version 2 
of \eqn{AD_{comb}}}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations used for each block of samples.}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, 
\emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}

\note{
This test is useful in analyzing treatment effects in randomized 
(incomplete) block experiments and in examining performance 
equivalence of several laboratories when presented with different 
test materials for comparison.
}
\seealso{
\code{\link{ad.test}}, \code{\link{ad.pval}}
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

## Run ad.test.combined.
set.seed(2627)
ad.test.combined(x1, x2, method = "simulated", Nsim = 1000) 
# or with same seed
# ad.test.combined(list(x1, x2), method = "simulated", Nsim = 1000)
# ad.test.combined(A~G|B,data=datx1x2,method="simulated",Nsim=1000)


}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
