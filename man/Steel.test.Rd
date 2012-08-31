\name{Steel.test}
\alias{Steel.test}
\title{
Steel's Mutiple Comparison Wilcoxon Tests
}
\description{
This function uses pairwise Wilcoxon tests, comparing a common control sample with each of
several treatment samples, in a multiple comparison fashion. The experiment wise
significance probabity is calculated, estimated, or approximated, when
testing
the hypothesis that all independent samples arise
from a common unspecified distribution, or that treatments have no effect when assigned
randomly to the given subjects.
}
\usage{
Steel.test(\dots, method = c("asymptotic", "simulated", "exact"),
	alternative = c("greater","less","two-sided"),
	continuity.corr = TRUE, dist = FALSE, Nsim = 10000)
}
\arguments{
  	\item{\dots}{
		Either several sample vectors, say 
		\eqn{x_1, \ldots , x_k}{x.1, \ldots , x.k}, 
		with \eqn{x_i}{x.i} containing \eqn{n_i}{n.i} sample values
		(\eqn{n_i > 4}{n.i > 4} is recommended for reasonable asymptotic 
		\eqn{P}-value calculation.) The first sample vector is taken as
		the control sample, while the others are samples under \eqn{k-1 \ge 1} different 
		treatments. 

		or a list of such sample vectors. 
	}
    \item{method}{= \code{c("asymptotic","simulated","exact")}, where

		\code{"asymptotic"} uses only an asymptotic normal approximation 
		to approximate the \eqn{P}-value, This calculation is always done.

		\code{"simulated"} uses \code{Nsim} simulated standardized
		Steel statistics based on random splits of the 
		pooled samples into samples of sizes 
		\eqn{n_1, \ldots, n_k}{n.1, \ldots, n.k}, to estimate the \eqn{P}-value.

		\code{"exact"} uses full enumeration of all sample splits with resulting 
		standardized Steel statistics to obtain the exact \eqn{P}-value. 
     		It is used only when \code{Nsim} is
      		at least as large as the number of full enumerations. Otherwise, \code{method}
		reverts to \code{"simulated"} using the given \code{Nsim}.
    }
    \item{alternative}{= \code{c("greater","less","two-sided")}, where for \code{"greater"} the 
		maximum of the pairwise standardized Wilcoxon test statistics is used and 
                a large maximum value is judged significant.
		For \code{"less"} the minimum of the pairwise standardized Wilcoxon test 
		statistics is used and a low minimum value is judged significant.
		For \code{"two-sided"} the maximum of the absolute pairwise standardized Wilcoxon test 
		statistics is used and a large maximum value is judged significant.
		
    }
    \item{continuity.corr}{\code{= TRUE} (default) or \code{FALSE}. This indicates whether a 
  		continuity correction should take place in the normal approximation, and it 
		takes place only when there are no ties in the pooled samples.
    }
    \item{dist}{\code{= FALSE} (default) or \code{TRUE}. If \code{TRUE}, the 
		simulated or fully enumerated null distribution vector \code{null.dist} 
		is returned for the Steel test statistic, as chosen via \code{alternative}.
    }
    \item{Nsim}{\code{= 10000} (default), number of simulation sample splits to use.	
		It is only used when \code{method = "simulated"},
		or when \code{method = "exact"} reverts to \code{method = "simulated"}, 
		as previously explained.
    }

}
\details{
The Steel criterion uses the Wilcoxon test statistic in the pairwise comparisons of the 
common control sample with each of the treatment samples. These statistics are used in 
standardized form, using the means and standard deviations as they apply conditionally
given the tie pattern in the pooled data, see Scholz (2012). This conditional treatment allows for 
correct usage in the presence of ties and is appropriate either when the samples are independent
and come from the same distribution (continuous or not) or when treatments are assigned 
randomly among the total of \code{N} subjects. However, in the case of ties the significance probability
has to be viewed conditionally given the tie pattern.

The Steel statistic is used to test the hypothesis that the samples all come 
from the same but unspecified distribution function \eqn{F(x)}, or, under random treatment 
assigment, that the treatments have no effect. The significance probability is the probability
of obtaining test results as extreme or more extreme than the observed test statistic,
when testing for the possibility of a treatment effect under any of the treatments.

For small sample sizes exact (conditional) null distribution
calculations are possible (with or without ties), based on a recursively extended
version of Algorithm C (Chase's sequence) in Knuth (2011), which allows the 
enumeration of all possible splits of the pooled data into samples of
sizes of \eqn{n_1, \ldots, n_k}{n.1, \ldots, n.k}, as appropriate under treatment randomization. This 
is done in C, as is the simulation of such splits.

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.

}
\value{
A list of class \code{kSamples} with components 
\item{test.name}{\code{"Steel"}}
\item{k}{number of samples being compared, including the control sample as the first one}
\item{ns}{vector \eqn{(n_1,\ldots,n_k)}{(n.1, \ldots, n.k)} of the \eqn{k} sample sizes}
\item{N}{size of the pooled sample \eqn{= n_1+\ldots+n_k}{= n.1+\ldots+n.k}}
\item{n.ties}{number of ties in the pooled sample}
\item{st}{2 (or 3) vector containing the observed standardized Steel statistic, 
its asymptotic \eqn{P}-value, 
(its simulated or exact \eqn{P}-value)}
\item{warning}{logical indicator, \code{warning = TRUE} when at least one 
\eqn{n_i < 5}{n.i < 5}}
\item{null.dist}{simulated or enumerated null distribution vector
of the test statistic. It is \code{NULL} when \code{dist = FALSE} or when
\code{method = "asymptotic"}.}
\item{method}{the \code{method} used.}
\item{Nsim}{the number of simulations used.}
}

\references{
Knuth, D.E. (2011), \emph{The Art of Computer Programming, Volume 4A 
Combinatorial Algorithms Part 1}, Addison-Wesley

Lehmann, E.L. (2006),
\emph{Nonparametrics, Statistical Methods Based on Ranks, Revised First Edition},
Springer Verlag.

Scholz, F.W. (2012), "On Steel's Test with Ties", submitted.
}
\section{warning}{\code{method = "exact"} should only be used with caution.
Computation time is proportional to the number of enumerations. 
Experiment with \code{\link{system.time}} and trial values for
\code{Nsim} to get a sense of the required computing time.
In most cases
\code{dist = TRUE} should not be used, i.e.,  
when the returned distribution objects 
become too large for R's work space.}


\examples{
z1 <- c(103, 111, 136, 106, 122, 114)
z2 <- c(119, 100,  97,  89, 112,  86)
z3 <- c( 89, 132,  86, 114, 114, 125)
z4 <- c( 92, 114,  86, 119, 131,  94)
st.out <- Steel.test(list(z1,z2,z3,z4),method="simulated",   
 	alternative="less",continuity.corr=TRUE,Nsim=10000)
# or
# st.out <- Steel.test(z1,z2,z3,z4,method="simulated",   
#	alternative="less",continuity.corr=TRUE,Nsim=10000)
# names(out.st)
}

\keyword{nonparametric}
\keyword{htest}
\keyword{design}
