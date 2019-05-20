\name{ad.test}
\alias{ad.test}
\title{
Anderson-Darling k-Sample Test
}
\description{
This function uses the Anderson-Darling criterion to test
the hypothesis that \eqn{k} independent samples with sample sizes  
\eqn{n_1,\ldots, n_k} arose
from a common unspecified distribution function \eqn{F(x)} and testing is 
done conditionally given the observed tie pattern. Thus this is a permutation test.
Both versions of the \eqn{AD} statistic are computed.
}
\usage{
ad.test(\dots, data = NULL, method = c("asymptotic", "simulated", "exact"),
	dist = FALSE, Nsim = 10000)
}
\arguments{
	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots, x_k}, 
		with \eqn{x_i} containing \eqn{n_i} sample values.
		\eqn{n_i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation. The pooled sample size is denoted 
		by \eqn{N=n_1+\ldots+n_k},

		or a list of such sample vectors,

		or a formula y ~ g, where y contains the pooled sample values 
		and g is a factor (of same length as y) with levels identifying 
		the samples to which the elements of y belong. 
	}

    	\item{data}{= an optional data frame providing the variables in formula y ~ g.
		}

    	\item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic \eqn{P}-value approximation, reasonable 
		for P in [.00001, .99999] when all \eqn{n_i > 4}. 
		Linear extrapolation via \eqn{\log(P/(1-P))}{log(P/(1-P))}
		is used outside [.00001, .99999]. This calculation is always done.
		See \code{\link{ad.pval}} for details. 
		The adequacy of the asymptotic \eqn{P}-value calculation
	 	may be checked using \code{\link{pp.kSamples}}.

		\code{"simulated"} uses \code{Nsim} simulated \eqn{AD} statistics, based on random 
		splits of the pooled samples into samples of sizes 
		\eqn{n_1, \ldots, n_k}, to estimate the exact conditional \eqn{P}-value.

		\code{"exact"} uses full enumeration of all sample splits with 
		resulting \eqn{AD} statistics to obtain the exact conditional \eqn{P}-values. 
		It is used only when \code{Nsim} is at least as large as the number
		\deqn{ncomb = \frac{N!}{n_1!\ldots n_k!}}{N!/(n_1!\ldots n_k!)}
		of full enumerations. Otherwise, \code{method}
		reverts to \code{"simulated"} using the given \code{Nsim}. It also reverts
		to \code{"simulated"} when \eqn{ncomb > 1e8} and \code{dist = TRUE}.
        }
	\item{dist}{\code{= FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated distribution vectors \code{null.dist1} and
		\code{null.dist2} are returned for the respective test statistic versions.
		Otherwise, \code{NULL} is returned. When \code{dist = TRUE} then 
		\code{Nsim <- min(Nsim, 1e8)}, to limit object size.
	}
	\item{Nsim}{\code{= 10000} (default), number of simulation sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method =}
		\code{ "simulated"}, as previously explained.
	}

}
\details{
If \eqn{AD} is the Anderson-Darling criterion for the \eqn{k} samples, 
its standardized test statistic is \eqn{T.AD = (AD - \mu)/\sigma}, with 
\eqn{\mu = k-1} and
\eqn{\sigma} representing mean and standard deviation of \eqn{AD}. This statistic 
is used to test the hypothesis that the samples all come 
from the same but unspecified continuous distribution function \eqn{F(x)}.


According to the reference article, two versions
of the \eqn{AD} test statistic are provided.
The above mean and standard deviation are strictly
valid only for version 1 in the 
continuous distribution case.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

The continuity assumption can be dispensed with, if we deal with 
independent random samples, or if randomization was used in allocating
subjects to samples or treatments, and if we view
the simulated or exact \eqn{P}-values conditionally, given the tie pattern
in the pooled samples. Of course, under such randomization any conclusions 
are valid only with respect to the group of subjects that were randomly allocated
to their respective samples.
The asymptotic \eqn{P}-value calculation assumes distribution continuity. No adjustment 
for lack thereof is known at this point. For details on the asymptotic 
\eqn{P}-value calculation see \code{\link{ad.pval}}.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Anderson-Darling"}}
\item{k}{number of samples being compared}
\item{ns}{vector of the \eqn{k} sample sizes \eqn{(n_1,\ldots,n_k)}}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}}
\item{n.ties}{number of ties in the pooled samples}
\item{sig}{standard deviations \eqn{\sigma} of version 1 of \eqn{AD} under the continuity assumption}
\item{ad}{2 x 3 (2 x 4) matrix containing \eqn{AD, T.AD}, asymptotic \eqn{P}-value, 
(simulated or exact \eqn{P}-value), for each version of the standardized test statistic \eqn{T},
version 1 in row 1, version 2 in row 2.}
\item{warning}{logical indicator, warning = TRUE when at least one 
\eqn{n_i < 5}{n.i < 5}}
\item{null.dist1}{simulated or enumerated null distribution of version 1 
of the test statistic, given as vector of all generated \eqn{AD} statistics.}
\item{null.dist2}{simulated or enumerated null distribution of version 2 
of the test statistic, given as vector of all generated \eqn{AD} statistics.}
\item{method}{The \code{method} used.}
\item{Nsim}{The number of simulations.}
}

\section{warning }{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. In most cases
\code{dist = TRUE} should not be used, i.e.,  
when the returned distribution vectors \code{null.dist1} and \code{null.dist2}
become too large for the R work space. These vectors are limited in length by 1e8. 
}


\note{
For small sample sizes and small \eqn{k} exact null distribution
calculations are possible (with or without ties), based on a recursively extended
version of Algorithm C (Chase's sequence) in Knuth (2011), Ch. 7.2.1.3, which allows the 
enumeration of all possible splits of the pooled data into samples of
sizes of \eqn{n_1, \ldots, n_k}, as appropriate under treatment randomization. The
enumeration and simulation are both done in C.
}

\note{
It has recently come to our attention that the Anderson-Darling test, originally 
proposed by Pettitt (1976) in the 2-sample case and generalized to k samples by
Scholz and Stephens, has a close relative created by Baumgartner et al (1998) in
the 2 sample case and populatized by Neuhaeuser (2012) with at least 6 papers 
among his cited references and generalized by Murakami (2006) to k samples.
}

\references{
Baumgartner, W., Weiss, P. and Schindler, H. (1998), A nonparametric test for the 
general two-sample problem, \emph{Bionetrics}, \bold{54}, 1129-1135.

Knuth, D.E. (2011), \emph{The Art of Computer Programming, Volume 4A 
Combinatorial Algorithms Part 1}, Addison-Wesley


Neuhaeuser, M. (2012), \emph{Nonparametric Statistical Tests, A Computational Approach},
CRC Press.

Murakami, H. (2006), A k-sample rank test based on modified Baumgartner statistic and
it power comparison, 
\emph{Jpn. Soc. Comp. Statist.}, \bold{19}, 1-13.

Murakami, H. (2012), Modified Baumgartner statistic for the two-sample and multisample 
problems: a numerical comparison.
\emph{J. of Statistical Comput. and Simul.}, \bold{82:5}, 711-728.

Pettitt, A.N. (1976), A two-sample Anderson_Darling rank statistic, \emph{Biometrika},
\bold{63}, 161-168.

Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, 
\emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}

\seealso{
\code{\link{ad.test.combined}}, \code{\link{ad.pval}}
}
\examples{
u1 <- c(1.0066, -0.9587,  0.3462, -0.2653, -1.3872)
u2 <- c(0.1005, 0.2252, 0.4810, 0.6992, 1.9289)
u3 <- c(-0.7019, -0.4083, -0.9936, -0.5439, -0.3921)
y <- c(u1, u2, u3)
g <- as.factor(c(rep(1, 5), rep(2, 5), rep(3, 5)))
set.seed(2627)
ad.test(u1, u2, u3, method = "exact", dist = FALSE, Nsim = 1000)
# or with same seed
# ad.test(list(u1, u2, u3), method = "exact", dist = FALSE, Nsim = 1000)
# or with same seed
# ad.test(y ~ g, method = "exact", dist = FALSE, Nsim = 1000)
}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
