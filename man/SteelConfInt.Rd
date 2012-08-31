\name{SteelConfInt}
\alias{SteelConfInt}
\title{
Simultaneous Confidence Bounds Based on Steel's Mutiple Comparison Wilcoxon Tests
}
\description{
This function inverts pairwise Wilcoxon tests, comparing a common control sample with each of
several treatment samples, in a multiple comparison fashion, to provide simultaneous confidence bounds
for the respective shift parameters by which the sampled continuous populations may differ
from each other. The joint coverage probability is calculated, estimated, or approximated, see
details.
}
\usage{
SteelConfInt(\dots, conf.level = 0.95, 
	alternative = c("less", "greater", "two.sided"), 
     	method = c("asymptotic", "exact", "simulated"),
	Nsim = 10000, cc = 0.5)
}
\arguments{
  	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots , x_k}{x.1, \ldots , x.k}, 
		with \eqn{x_i}{x.i} containing \eqn{n_i}{n.i} sample values
		(\eqn{n_i > 4}{n.i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation.) The first sample vector is taken as
		the control sample, while the others are samples under \eqn{s=k-1} different 
		treatments. 

		or a list of such sample vectors. 
	}
    \item{conf.level}{\code{= 0.95} (default) the target joint confidence level for all bounds/intervals.

		\code{0 < conf.level < 1}.
    }
    \item{alternative}{= \code{c("less", "greater", "two.sided")}, where \code{"less"} results in
		simultaneous upper confidence bounds for all shift parameters and any negative values should 
		lead to the rejection of the null hypothesis of the corresponding shift differences 
		being zero in favor of at least one being less than zero.

		\code{"greater"} results in simultaneous lower confidence bounds for all shift parameters 
		and any positive values should lead to the rejection of the same null hypothesis in favor of 
		at least one being greater than zero.

		\code{"two.sided"} results in simultaneous confidence intervals for all shift parameters 
		and any interval not straddling zero should lead to the rejection of the same null hypothesis 
		in favor of at least one being different from zero.
    }
    \item{method}{= \code{c("asymptotic", "exact", "simulated")}, where

		\code{"asymptotic"} uses only an asymptotic normal approximation 
		to approximate the achieved coverage probability. 
		This calculation is always done.

		\code{"exact"} uses full enumeration of all sample splits to obtain the exact achieved 
		coverage probability (see Details).
     		It is used only when \code{Nsim} is at least as large as the number of full enumerations. 
		Otherwise, \code{method} reverts to \code{"simulated"} using the given \code{Nsim}.

		\code{"simulated"} uses \code{Nsim} simulated random splits of the 
		pooled samples into samples of sizes \eqn{n_1, \ldots, n_k}{n.1, \ldots, n.k}, 
		to estimate the achieved coverage probability.


    }
    \item{Nsim}{\code{= 10000} (default), number of simulation sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method = "simulated"}, 
		as previously explained.
    }
    \item{cc}{\code{= 0.5} (default), should be 0 or 0.5. The latter uses a continuity correction in the 
		normal approximation, the former does not.
    }

}
\details{
Let \eqn{W_{0i}, i=1,\ldots,s}{W.0i, i=1,\ldots,s} denote the respective Wilcoxon statistics comparing the common control sample
with each of the \eqn{s} treatment samples. Furthermore, let \eqn{D_{ij}}{D.ij} denote the \eqn{j}-th largest of 
all paired differences between the observations in the \eqn{i}-th treatment sample and those of the control
sample. By simple extension of results in Lehmann (2006) the following equation holds, relating the null distribution of the 
Wilcoxon statistics and the coverage probabilities of the \eqn{D_{ij_i}}{D.ij.i} for any set of
\eqn{j_1,\ldots,j_s}{j11,\ldots,j.s}
\deqn{P_\Delta(\Delta_i \le D_{ij_i}, i=1,\ldots,s)=P_0(W_{0i}\le j_i -1, i=1,\ldots,s)}{P.\Delta(\Delta.i \le D.ij.i, i=1,\ldots,s)=P.0(W.0i\le j.i -1, i=1,\ldots,s)}
and
\deqn{P_\Delta(\Delta_i \ge D_{ij_i}, i=1,\ldots,s)=P_0(W_{0i}\le n_0 n_i -j_i, i=1,\ldots,s)}{P.\Delta(\Delta.i \ge D.ij_i, i=1,\ldots,s)=P.0(W.0i\le n.0 n.i -j.i, i=1,\ldots,s)}
where \eqn{P_\Delta}{P.\Delta} refers to the distribution under \eqn{\Delta=(\Delta_1,\ldots,\Delta_s)}{\Delta=(\Delta.1,\ldots,\Delta.s)}
and \eqn{P_0}{P.0} refers to the joint null distribution of the \eqn{W_{0i}}{W.0i} when all sampled 
distributions are the same and continuous. There are \eqn{s} indices \eqn{j_i}{j.i} that can be manipulated
to affect the achieved confidence level. To limit the computational complexity
standardized versions of the \eqn{W_{0i}}{W.0i}, i.e.,  \eqn{(W_{0i}-\mu_i)/\tau_i}{(W.0i-\mu.i)/\tau.i} with 
\eqn{\mu_i}{\mu.i} and \eqn{\tau_i}{\tau.i} representing mean and standard deviation of \eqn{W_{0i}}{W.0i},
are used to choose a common value for \eqn{(j_i -1-\mu_i)/\tau_i}{(j.i -1-\mu.i)/\tau.i}  from the multivariate normal approximation 
for the \eqn{W_{0i}} (see Miller (1981)), and reduce that 
to integer values for \eqn{j_i}{j.i}, rounding up, rounding down, and rounding to the nearest integer. These
integers  \eqn{j_i}{j.i} are then used in approximating the actual probabilities
\eqn{P_0(W_{0i}\le j_i -1, i=1,\ldots,s)}{P.0(W.0i\le j.i -1, i=1,\ldots,s)}, and from these three coverage probabilities 
the one that is closest to the nominal confidence level \eqn{\gamma} and \eqn{\ge \gamma}
and also also the one that is closest without the restriction \eqn{\ge \gamma} are chosen.

When \code{method = "exact"} or \code{= "simulated"} is specified, the same process
is used, using either the fully enumerated exact distribution of \eqn{W_{0i}, i=1,\ldots,s}{W.0i, i=1,\ldots,s} (based on a recursive 
version of Chase's sequence as presented in Knuth (2011)) for all sample splits,
or the simulated distribution of \eqn{W_{0i}, i=1,\ldots,s}{W.0i, i=1,\ldots,s}. However, since these distribtions are discrete
the starting point before rounding up is the smallest quantile such that the proportion of distribution values less 
or equal to it is at least \eqn{\gamma}. The starting point before rounding down is the highest quantile such that 
the proportion of distribution values less 
or equal to it is at most \eqn{\gamma}. The third option of rounding to the closest integer is performed using 
the average of the first two.

Confidence intervals are constructed by using upper and lower confidence bounds, each with
same confidence level of \eqn{(1+\gamma)/2}.

When the orginal sample data appear to be rounded, and especially when there are ties,
one should widen the computed intervals or bounds by the rounding \eqn{\epsilon}, as illustrated
in Lehmann (2006). For example, when all sample values appear to end in one of \eqn{.0, .2, .4, .6, .8},
the rounding \eqn{\epsilon} would be \eqn{.2}. Ultimately, this is a judgment call for the user. Such widening
of intervals will make the actually achieved confidence level \eqn{\ge} the stated achieved level. 
}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Steel.bounds"}}
\item{n0}{the control sample size \eqn{= n_1}{n.1} }
\item{ns}{vector \eqn{(n_2,\ldots,n_k)}{(n.2, \ldots, n.k)} of the \eqn{s=k-1} treatment sample sizes}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}{= n.1+\ldots+n.k}}
\item{n.ties}{number of ties in the pooled sample}
\item{bounds}{a list of data frames. When \code{method = "asymptotic"} is specified, the list
consists of two data frames named \code{conservative.bounds.asymptotic}
and \code{closest.bounds.asymptotic}. Each data frame consists of \eqn{s} rows
corresponding to the \eqn{s} shift parameters \eqn{\Delta_i}{\Delta.i} and three columns,
the first column giving the lower bound, the second column the upper bound, while
the first row of the third column states the computed confidence level by asymptotic 
approximation, applying jointly to all \eqn{s} sets of bounds. For one-sided bounds
the corresponding other bound is set to \code{Inf} or \code{-Inf}, respectively.

In case of \code{conservative.bounds.asymptotic} the achieved asymptotic confidence level is 
targeted to be \eqn{\ge} \code{conf.level}, but closest to it among three possible choices (see details).

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
\item{j.LU}{an \eqn{s} by 4 matrix giving the indices \eqn{j_i}{j.i} used for computing the bounds
\eqn{D_{ij_i}}{D.ij.i} for \eqn{\Delta_i, i=1,\ldots, s}{\Delta.i, i=1,\ldots, s}.
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
SteelConfInt(list(z1,z2,z3,z4),conf.level=0.95,alternative="two.sided",
	method="simulated",Nsim=10000)
# or with same results
# SteelConfInt(z1,z2,z3,z4,conf.level=0.95,alternative="two.sided",
#	method="simulated",Nsim=10000)
}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
