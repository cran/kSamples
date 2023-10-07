\name{kSamples-package}
\alias{kSamples-package}
\alias{kSamples}
\docType{package}
\title{
The Package kSamples Contains Several Nonparametric K-Sample Tests and 
their Combinations over Blocks
}
\description{
The k-sample Anderson-Darling, Kruskal-Wallis, normal score and van der Waerden score tests
are used to test the hypothesis that k samples of 
sizes \eqn{n_1, \ldots, n_k}
come from a common continuous distribution \eqn{F} that is otherwise unspecified. They are rank tests.
Average rank scores are used in case of ties. 
While \code{\link{ad.test}} is consistent against all alternatives, \code{\link{qn.test}}
tends to be sensitive mainly to shifts between samples.
The combined versions of these tests,
\code{\link{ad.test.combined}} and
\code{\link{qn.test.combined}}, are 
used to simultaneously test such hypotheses across several blocks of samples.
The hypothesized common distributions and the number k of samples for each block of samples
may vary from block to block. 

The Jonckheere-Terpstra test addresses the same hypothesis as above but is sensitive to increasing
alternatives (stochastic ordering).

Also treated is the analysis of 2 x t contingency tables
using the Kruskal-Wallis criterion and its extension to blocks.
 
Steel's simultaneous comparison test of a common control sample with \eqn{s=k-1} treatment samples
using pairwise Wilcoxon tests for each control/treatment pair is provided, and also 
the simultaneous confidence bounds of treatment shift effects resulting from the inversion of these tests
when sampling from continuous populations.

Distributional aspects are handled asymptotically in all cases, and by choice also
via simulation or exact enumeration.
While simulation is always an option, exact calculations
are only possible for small sample sizes and only when few samples are involved. These exact
calculations can be done with or without ties in the pooled samples, based on a recursively extended
version of Algorithm C (Chase's sequence) in Knuth (2011), which allows the 
enumeration of all possible splits of the pooled data into samples of
sizes of \eqn{n_1, \ldots, n_k}, as appropriate under treatment randomization
or random sampling, when viewing tests conditionally given the observed tie pattern.
}


\author{
Fritz Scholz and Angie Zhu

Maintainer: Fritz Scholz <fscholz@u.washington.edu>
}
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

Lehmann, E.L. (2006), \emph{Nonparametrics, Statistical Methods Based on Ranks},
Revised First Edition, 
Springer, New York.

Scholz, F.W. (2023), "On Steel's Test with Ties", \url{https://arxiv.org/abs/2308.05873}

Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, \emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}


\keyword{nonparametric}
\keyword{htest}
\keyword{design}
