\name{JT.dist}
\alias{djt}
\alias{pjt}
\alias{qjt}
\title{
Null Distribution of the 
Jonckheere-Terpstra k-Sample Test Statistic
}
\description{
The
Jonckheere-Terpstra k-sample test statistic JT is defined
as \eqn{JT = \sum_{i<j} W_{ij}} where
\eqn{W_{ij}}{W_{ij}} is the Mann-Whitney statistic comparing 
samples \eqn{i} and \eqn{j}, indexed in the order 
of the stipulated increasing alternative. 
It is assumed that there are no ties
in the pooled samples.

This function uses Harding's algorithm as far as computations
are possible without becoming unstable.
}
\usage{
djt(x, nn)

pjt(x, nn)

qjt(p, nn)
}
\arguments{
  	\item{x}{a numeric vector, typically integers}
		
	\item{nn}{a vector of integers, representing the sample sizes
		in the order stipulated by the alternative}

	\item{p}{a vector of probabilities}

}
\details{While Harding's algorithm is mathematically correct, 
it is problematic in its computing implementation. 
The counts become very large and normalizing 
them by combinatorials leads to significance loss. 
When that happens the functions return an error message: 
can't compute due to numerical instability.
This tends to happen when the total number of sample values 
becomes too large.
That depends also on the way the sample sizes are allocated. 
}
\value{
For \code{djt} it is a vector 
\eqn{p = (p_1,\ldots,p_n)}
giving the values of 
\eqn{p_i = P(JT = x_i)}, where \code{n} is the length
of the input \code{x}.

For \code{pjt} it is a vector 
\eqn{P = (P_1,\ldots,P_n)}
giving the values of 
\eqn{P_i = P(JT \leq x_i)}{P_i = P(JT <= x_i)}.

For \code{qjt} is a vecto r \eqn{x = (x_1,\ldots,x_n)},where \eqn{x_i}
is the smallest \eqn{x}
such that \eqn{P(JT \leq x) \geq p_i}{P(JT <= x) >= p_i}.


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
djt(c(-1.5,1.2,3), 2:4)
pjt(c(2,3.4,7), 3:5)
qjt(c(0,.2,.5), 2:4)
}
\keyword{nonparametric}
\keyword{htest}
\keyword{design}
