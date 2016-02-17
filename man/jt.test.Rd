\name{jt.test}
\alias{jt.test}

\title{
Jonckheere-Terpstra k-Sample Test for Increasing Alternatives
}
\description{
The
Jonckheere-Terpstra k-sample test statistic JT is defined
as \eqn{JT = \sum_{i<j} W_{ij}} where
\eqn{W_{ij}} is the Mann-Whitney statistic comparing 
samples \eqn{i} and \eqn{j}, indexed in the order 
of the stipulated increasing alternative. 
There may be ties in the pooled samples.
}
\usage{
jt.test(\dots, data = NULL, method=c("asymptotic","simulated","exact"),
		dist = FALSE, Nsim = 10000) 
}
\arguments{
	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots, x_k}, 
		with \eqn{x_i} containing \eqn{n_i} sample values.
		\eqn{n_i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation. The pooled sample size is denoted 
		by \eqn{N=n_1+\ldots+n_k}. The order of samples should be
		as stipulated under the alternative

		or a list of such sample vectors,

		or a formula y ~ g, where y contains the pooled sample values 
		and g (same length as y) is a factor with levels identifying 
		the samples to which the elements of y belong, the factor levels
		reflecting the order under the stipulated alternative,

	}

    	\item{data}{= an optional data frame providing the variables in formula y ~ g.
		}
	
    	\item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic normal \eqn{P}-value approximation.

		\code{"simulated"} uses \code{Nsim} simulated \eqn{JT} statistics based on random splits of the 
		pooled samples into samples of sizes 
		\eqn{n_1, \ldots, n_k}, to estimate the \eqn{P}-value.

		\code{"exact"} uses full enumeration of all sample splits with 
		resulting \eqn{JT} statistics to obtain the exact \eqn{P}-value. 
		It is used only when \code{Nsim} is at least as large as the number
		\deqn{ncomb = \frac{N!}{n_1!\ldots n_k!}}{N!/(n_1!\ldots n_k!)}
		of full enumerations. Otherwise, \code{method}
		reverts to \code{"simulated"} using the given \code{Nsim}. It also reverts
		to \code{"simulated"} when \eqn{ncomb > 1e8} and \code{dist = TRUE}.
    	}
	\item{dist}{\code{= FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated distribution vector \code{null.dist}
		is returned for the JT test statistic. Otherwise, \code{NULL} is returned. 
		When \code{dist = TRUE} then \code{Nsim <- min(Nsim, 1e8)}, 
		to limit object size.
	}
	\item{Nsim}{\code{= 10000} (default), number of simulation sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method =}
		\code{ "simulated"}, as previously explained.
	}

}
\details{
The JT statistic 
is used to test the hypothesis that the samples all come 
from the same but unspecified continuous distribution function \eqn{F(x)}.
It is specifically aimed at alternatives where the sampled distributions
are stochastically increasing.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

The continuity assumption can be dispensed with, if we deal with 
independent random samples, or if randomization was used in allocating
subjects to samples or treatments, and if we view
the simulated or exact \eqn{P}-values conditionally, given the tie pattern
in the pooled samples. Of course, under such randomization any conclusions 
are valid only with respect to the group of subjects that were randomly allocated
to their respective samples.
The asymptotic \eqn{P}-value calculation is valid provided all sample sizes become large.
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Jonckheere-Terpstra"}}
\item{k}{number of samples being compared}
\item{ns}{vector \eqn{(n_1,\ldots,n_k)} of the \eqn{k} sample sizes}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}}
\item{n.ties}{number of ties in the pooled sample}
\item{qn}{4 (or 5) vector containing the observed \eqn{JT}, its mean and standard deviation
and its asymptotic \eqn{P}-value, 
(and its simulated or exact \eqn{P}-value)}
\item{warning}{logical indicator, \code{warning = TRUE} when at least one 
\eqn{n_i < 5}{n.i < 5}}
\item{null.dist}{simulated or enumerated null distribution 
of the test statistic. It is \code{NULL} when \code{dist = FALSE} or when
\code{method = "asymptotic"}.}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations used.}

}


\references{
Harding, E.F. (1984), An Efficient, Minimal-storage Procedure for
Calculating the Mann-Whitney U, Generalized U and Similar Distributions,
\emph{Appl. Statist.} \bold{33} No. 1, 1-6.


Jonckheere, A.R. (1954), A Distribution Free \emph{k}-sample Test against Ordered 
Alternatives,
\emph{Biometrika}, \bold{41}, 133-145.

Lehmann, E.L. (2006),
\emph{Nonparametrics, Statistical Methods Based on Ranks, Revised First Edition},
Springer Verlag.

Terpstra, T.J. (1952), The Asymptotic Normality and Consistency of Kendall's Test against Trend, when Ties are Present in One Ranking,
\emph{Indagationes Math.} \bold{14}, 327-333.
}




\examples{
x1 <- c(1,2)
x2 <- c(1.5,2.1)
x3 <- c(1.9,3.1)
yy <- c(x1,x2,x3)
gg <- as.factor(c(1,1,2,2,3,3))
jt.test(x1, x2, x3,method="exact",Nsim=90)
# or 
# jt.test(list(x1, x2, x3), method = "exact", Nsim = 90)
# or
# jt.test(yy ~ gg, method = "exact", Nsim = 90)
}
\keyword{nonparametric}
\keyword{htest}
\keyword{design}
