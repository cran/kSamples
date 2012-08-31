\name{qn.test}
\alias{qn.test}
\title{
Rank Score k-Sample Tests
}
\description{
This function uses the \eqn{QN} criterion (Kruskal-Wallis, van der Waerden scores, normal scores) to test
the hypothesis that \eqn{k} independent samples arise
from a common unspecified distribution.
}
\usage{
qn.test(\dots, test = c("KW", "vdW", "NS"), 
	method = c("asymptotic", "simulated", "exact"),
	dist = FALSE, Nsim = 10000)
}
\arguments{
  	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots , x_k}{x.1, \ldots , x.k}, 
		with \eqn{x_i}{x.i} containing \eqn{n_i}{n.i} sample values, with
		\eqn{N=n_1+\ldots+n_k}{N=n.1+\ldots+n.k} denoting the pooled sample size.
		\eqn{n_i > 4}{n.i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation.

		or a list of such sample vectors. 
	}
	\item{test}{= \code{c("KW", "vdW", "NS")}, where

		\code{"KW"} uses scores \code{1:N} (Kruskal-Wallis test)

		\code{"vdW"} uses van der Waerden scores, \code{qnorm( (1:N) / (N+1) )}

		\code{"NS"} uses normal scores, i.e., expected standard normal order statistics,
		invoking function \code{normOrder} of \code{package SuppDists}
	}
    \item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic chi-square approximation 
		with \code{k-1} degrees of freedom to approximate the \eqn{P}-value. 
		This calculation is always done.

		\code{"simulated"} uses \code{Nsim} simulated \eqn{QN} statistics based on random splits of the 
		pooled samples into samples of sizes 
		\eqn{n_1, \ldots, n_k}{n.1, \ldots, n.k}, to estimate the \eqn{P}-value.

		\code{"exact"} uses full enumeration of all sample splits with resulting 
		\eqn{QN} statistics to obtain the exact \eqn{P}-value. 
     	It is used only when \code{Nsim} is
      	at least as large as the number of full enumerations. Otherwise, \code{method}
		reverts to \code{"simulated"} using the given \code{Nsim}.
    }
	\item{dist}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated null distribution vector \code{null.dist} 
		is returned for the \eqn{QN} test statistic.
	}
	\item{Nsim}{\code{= 10000} (default), number of simulation sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method = "simulated"}, 
		as previously explained.
	}

}
\details{
The \eqn{QN} criterion based on rank scores \eqn{v_1,\ldots,v_N}{v.1,\ldots,v.N} is
\deqn{QN=\frac{1}{s_v^2}\left(\sum_{i=1}^k \frac{(S_{iN}-n_i \bar{v}_{N})^2}{n_i}\right)}{
[ sum (S.iN - n.i * vbar.N)^2/n.i ] / s.v^2} 
where \eqn{S_{iN}}{S.iN} is the sum of rank scores for the \eqn{i}-th sample and 
\eqn{\bar{v}_N}{vbar.N} and 
\eqn{s_v^2}{s.v^2} are sample mean and sample variance (denominator \eqn{N-1})
of all scores.

The statistic \eqn{QN} is used to test the hypothesis that the samples all come 
from the same but unspecified continuous distribution function \eqn{F(x)}.
\eqn{QN} is always adjusted for ties by averaging the scores of tied observations.

Conditions for the asymptotic approximation (chi-square  with \eqn{k-1} degrees of freedom)  
can be found in Lehmann, E.L. (2006), Appendix Corollary 10, or in 
Hajek, Sidak, and Sen (1999), Ch. 6, problems 13 and 14.

For small sample sizes exact null distribution
calculations are possible (with or without ties), based on a recursively extended
version of Algorithm C (Chase's sequence) in Knuth (2011), which allows the 
enumeration of all possible splits of the pooled data into samples of
sizes of \eqn{n_1, \ldots, n_k}{n.1, \ldots, n.k}, as appropriate under treatment randomization. This 
is done in C, as is the simulation.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

The continuity assumption can be dispensed with, if we deal with 
independent random samples from any common distribution, 
or if randomization was used in allocating
subjects to samples or treatments, and if 
the asymptotic, simulated or exact \eqn{P}-values are viewed conditionally, given the tie pattern
in the pooled sample. Under such randomization any conclusions 
are valid only with respect to the subjects that were randomly allocated
to their respective treatment samples.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Kruskal-Wallis"}, \code{"van der Waerden scores"}, or

\code{"normal scores"}}
\item{k}{number of samples being compared}
\item{ns}{vector \eqn{(n_1,\ldots,n_k)}{(n.1, \ldots, n.k)} of the \eqn{k} sample sizes}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}{= n.1+\ldots+n.k}}
\item{n.ties}{number of ties in the pooled sample}
\item{qn}{2 (or 3) vector containing the observed \eqn{QN}, its asymptotic \eqn{P}-value, 
(its simulated or exact \eqn{P}-value)}
\item{warning}{logical indicator, \code{warning = TRUE} when at least one 
\eqn{n_i < 5}{n.i < 5}}
\item{null.dist}{simulated or enumerated null distribution 
of the test statistic. It is \code{NULL} when \code{dist = FALSE} or when
\code{method = "asymptotic"}.}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations used.}
}

\section{warning}{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. 
Experiment with \code{\link{system.time}} and trial values for
\code{Nsim} to get a sense of the required computing time.
In most cases
\code{dist = TRUE} should not be used, i.e.,  
when the returned distribution objects 
become too large for R's work space.}

\references{
Hajek, J., Sidak, Z., and Sen, P.K. (1999), \emph{Theory of Rank Tests (Second Edition)}, Academic Press.

Knuth, D.E. (2011), \emph{The Art of Computer Programming, Volume 4A 
Combinatorial Algorithms Part 1}, Addison-Wesley

Kruskal, W.H. (1952), A Nonparametric Test for the Several Sample Problem,
\emph{The Annals of Mathematical Statistics},
\bold{Vol 23, No. 4}, 525-540

Kruskal, W.H. and Wallis, W.A. (1952), Use of Ranks in One-Criterion Variance Analysis,
\emph{Journal of the American Statistical Association}, 
\bold{Vol 47, No. 260}, 583--621. 

Lehmann, E.L. (2006),
\emph{Nonparametrics, Statistical Methods Based on Ranks, Revised First Edition},
Springer Verlag.
}


\seealso{
\code{\link{qn.test.combined}}
}

\examples{
u1 <- c(1.0066, -0.9587,  0.3462, -0.2653, -1.3872)
u2 <- c(0.1005, 0.2252, 0.4810, 0.6992, 1.9289)
u3 <- c(-0.7019, -0.4083, -0.9936, -0.5439, -0.3921)
out.qn <- qn.test(u1,u2,u3,test="KW",method="simulated",dist=FALSE,Nsim=1000)
# or
# out.qn <- qn.test(list(u1,u2,u3),test="KW",method="simulated",
#		dist=FALSE,Nsim=1000)
# names(out.qn)
}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
