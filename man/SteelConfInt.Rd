\name{SteelConfInt}
\alias{SteelConfInt}
\title{
Simultaneous Confidence Bounds Based on Steel's Multiple Comparison Wilcoxon Tests
}
\description{
This function inverts pairwise Wilcoxon tests, comparing a common control sample with each of
several treatment samples to provide simultaneous confidence bounds
for the respective shift parameters by which the sampled treatment populations may differ
from the control population. It is assumed that all samples are independent 
and
that the sampled distributions are continuous to avoid ties.
The joint coverage probability for all bounds/intervals is calculated, estimated, or approximated, see
Details. For treatment of ties also see Details.
}
\usage{
SteelConfInt(\dots, data = NULL, conf.level = 0.95, 
	alternative = c("less", "greater", "two.sided"), 
     	method = c("asymptotic", "exact", "simulated"), Nsim = 10000)
}
\arguments{
	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots, x_k}, 
		with \eqn{x_i} containing \eqn{n_i} sample values.
		\eqn{n_i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation. The pooled sample size is denoted 
		by \eqn{N=n_1+\ldots+n_k}. The first vector serves as control sample,
		the others as treatment samples.

		or a list of such sample vectors.

		or a formula y ~ g, where y contains the pooled sample values 
		and g (same length as y) is a factor with levels identifying 
		the samples to which the elements of y belong. The lowest factor level
		corresponds to the control sample, the other levels to treatment samples. 
	}

    	\item{data}{= an optional data frame providing the variables in formula y ~ g.
		}

    \item{conf.level}{\code{= 0.95} (default) the target joint confidence level for all bounds/intervals.

		\code{0 < conf.level < 1}.
    }

    \item{alternative}{= \code{c("less", "greater", "two.sided")}, where \code{"less"} results in
		simultaneous upper confidence bounds for all shift parameters and any negative upper bound
		should lead to the rejection of the null hypothesis of all shift parameters 
		being zero or positive in favor of at least one being less than zero.

		\code{"greater"} results in simultaneous lower confidence bounds for all shift parameters 
		and any positive lower bound should lead to the rejection of the null hypothesis of all
		shift parameters being zero or negative in favor of at least one being greater 
		than zero.


		\code{"two.sided"} results in simultaneous confidence intervals for all shift parameters 
		and any interval not straddling zero should lead to the rejection of the null hypothesis 
		of all shift parameters being zero in favor of at least one being different from zero.
    }

    \item{method}{= \code{c("asymptotic", "exact", "simulated")}, where

		\code{"asymptotic"} uses only an asymptotic normal approximation 
		to approximate the achieved joint coverage probability. 
		This calculation is always done.

		\code{"exact"} uses full enumeration of all sample splits to obtain the exact achieved 
		joint coverage probability (see Details).
     		It is used only when \code{Nsim} is at least as large as the number of full enumerations. 
		Otherwise, \code{method} reverts to \code{"simulated"} using the given \code{Nsim}.

		\code{"simulated"} uses \code{Nsim} simulated random splits of the 
		pooled samples into samples of sizes \eqn{n_1, \ldots, n_k}, 
		to estimate the achieved joint coverage probability.


    }

    \item{Nsim}{\code{= 10000} (default), number of simulated sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method =}
		\code{"simulated"}, as previously explained.
    }

}
\details{The first sample is treated as control sample with sample size \eqn{n_1}. The remaining 
\eqn{s=k-1} samples are treatment samples.
Let \eqn{W_{1i}, i=2,\ldots,k} denote the respective Wilcoxon statistics comparing the common control sample (index 1)
with each of the \eqn{s} treatment samples (indexed by \eqn{i}). 
For each comparison of control and treatment \eqn{i}
sample
only the observations of the two samples involved are ranked.
By \eqn{W_i=W_{1i}-n_i(n_i+1)/2} we denote
the corresponding Mann-Whitney test statistic.
Furthermore, let \eqn{D_{i(j)}} denote the \eqn{j}-th ordered value (ascending order) of the \eqn{n_1n_i}
paired differences between the observations in treatment sample \eqn{i} and those of the control
sample. By simple extension of results in Lehmann (2006), pages 87 and 92, the following equations hold, 
relating the null distribution of the 
Mann-Whitney statistics and the joint coverage probabilities of the \eqn{D_{i(j_i)}} for any set of
\eqn{j_1,\ldots,j_s} with \eqn{1\le j_i \le n_1 n_i}.
\deqn{P_\Delta(\Delta_i \le D_{i(j_i)}, i=2,\ldots,k)=P_0(W_i\le j_i -1, i=2,\ldots,k)}
and
\deqn{P_\Delta(\Delta_i \ge D_{i(j_i)}, i=2,\ldots,s)=P_0(W_{i}\le n_1 n_i -j_i, i=2,\ldots,k)}
where \eqn{P_\Delta} refers to the distribution under \eqn{\Delta=(\Delta_2,\ldots,\Delta_k)}
and \eqn{P_0} refers to the joint null distribution of the \eqn{W_i} when all sampled 
distributions are the same and continuous. There are \eqn{k-1} indices \eqn{j_i} that can be manipulated
to affect the achieved confidence level. To limit the computational complexity
standardized versions of the \eqn{W_i}, i.e.,  \eqn{(W_i-\mu_i)/\tau_i} with 
\eqn{\mu_i} and \eqn{\tau_i} representing mean and standard deviation of \eqn{W_i},
are used to choose a common value for \eqn{(j_i -1-\mu_i)/\tau_i}  (satisfying the 
\eqn{\gamma} level) from the multivariate normal approximation 
for the \eqn{W_i} (see Miller (1981) and Scholz (2016)), and reduce that 
to integer values for \eqn{j_i}, rounding up, rounding down, and rounding to the nearest integer. These
integers  \eqn{j_i} are then used in approximating the actual joint probabilities
\eqn{P_0(W_i\le j_i -1, i=2,\ldots,k)}, and from these three coverage probabilities 
the one that is closest to the nominal confidence level \eqn{\gamma} and \eqn{\ge \gamma}
and also also the one that is closest without the restriction \eqn{\ge \gamma} are chosen.

When \code{method = "exact"} or \code{= "simulated"} is specified, the same process
is used, using either the fully enumerated exact distribution of \eqn{W_i, i=2,\ldots,k} (based on a recursive 
version of Chase's sequence as presented in Knuth (2011)) for all sample splits,
or the simulated distribution of \eqn{W_i, i=2,\ldots,k}. However, since these distributions are discrete
the starting point before rounding up is the smallest quantile such that the proportion of distribution values less 
or equal to it is at least \eqn{\gamma}. The starting point before rounding down is the highest quantile such that 
the proportion of distribution values less 
or equal to it is at most \eqn{\gamma}. The third option of rounding to the closest integer is performed using 
the average of the first two.

Confidence intervals are constructed by using upper and lower confidence bounds, each with
same confidence level of \eqn{(1+\gamma)/2}.

When the original sample data appear to be rounded, and especially when there are ties,
one should widen the computed intervals or bounds by the rounding \eqn{\epsilon}, as illustrated
in Lehmann (2006), pages 85 and 94. For example, when all sample values appear to end in one of \eqn{.0, .2, .4, .6, .8},
the rounding \eqn{\epsilon} would be \eqn{.2}. Ultimately, this is a judgment call for the user. Such widening
of intervals will make the actually achieved confidence level \eqn{\ge} the stated achieved level. 
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Steel.bounds"}}
\item{n1}{the control sample size \eqn{= n_1}}
\item{ns}{vector \eqn{(n_2,\ldots,n_k)} of the \eqn{s=k-1} treatment sample sizes}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}}
\item{n.ties}{number of ties in the pooled sample}
\item{bounds}{a list of data frames. When \code{method = "asymptotic"} is specified, the list
consists of two data frames named \code{conservative.bounds.asymptotic}
and \code{closest.bounds.asymptotic}. Each data frame consists of \eqn{s} rows
corresponding to the \eqn{s} shift parameters \eqn{\Delta_i} and three columns,
the first column giving the lower bound, the second column the upper bound, while
the first row of the third column states the computed confidence level by asymptotic 
approximation, applying jointly to all \eqn{s} sets of bounds. For one-sided bounds
the corresponding other bound is set to \code{Inf} or \code{-Inf}, respectively.

In case of \code{conservative.bounds.asymptotic} the achieved asymptotic confidence level is 
targeted to be \eqn{\ge} \code{conf.level}, but closest to it among three possible choices (see Details).

In the case of \code{closest.bounds.asymptotic} the achieved level is targeted to
be closest to \code{conf.level}.

When \code{method = "exact"} or \code{method = "simulated"}
is specified the list consists in addition of two further data frames named either

\code{conservative.bounds.exact} and 
\code{closest.bounds.exact} or

\code{conservative.bounds.simulated} and 
\code{closest.bounds.simulated}. 

In either case the structure and meaning
of these data frames parallels that of the \code{"asymptotic"} case. 
}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations used.}
\item{j.LU}{an \eqn{s} by 4 matrix giving the indices \eqn{j_i} used for computing the bounds
\eqn{D_{i(j_i)}} for \eqn{\Delta_i, i=1,\ldots, s}.
}
}

\references{
Knuth, D.E. (2011), \emph{The Art of Computer Programming, Volume 4A 
Combinatorial Algorithms Part 1}, Addison-Wesley

Lehmann, E.L. (2006),
\emph{Nonparametrics, Statistical Methods Based on Ranks, Revised First Edition},
Springer Verlag.

Miller, Rupert G., Jr. (1981), \emph{Simultaneous Statistical Inference, Second Edition},
Springer Verlag, New York.

Scholz, F.W. (2016), "On Steel's Test with Ties", submitted to \emph{Journal of Nonparametric Statistics}.
}
\section{warning}{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. 
Experiment with \code{\link{system.time}} and trial values for
\code{Nsim} to get a sense of the required computing time.}


\examples{
z1 <- c(103, 111, 136, 106, 122, 114)
z2 <- c(119, 100,  97,  89, 112,  86)
z3 <- c( 89, 132,  86, 114, 114, 125)
z4 <- c( 92, 114,  86, 119, 131,  94)
set.seed(2627)
SteelConfInt(list(z1,z2,z3,z4),conf.level=0.95,alternative="two.sided", 
   method="simulated",Nsim=10000)
# or with same seed
# SteelConfInt(z1,z2,z3,z4,conf.level=0.95,alternative="two.sided", 
#   method="simulated",Nsim=10000)
}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
