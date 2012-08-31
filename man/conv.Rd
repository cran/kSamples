\name{conv}
\alias{conv}
\title{
Convolution of Two Discrete Distributions
}
\description{
This function convolutes two discrete distribution, each given by strictly increasing support vectors
and corresponding probability vectors.
}
\usage{
	conv(x1,p1,x2,p2)
}
\arguments{
  	\item{x1}{
		support vector of the first distribution, with strictly increasing elements.
	}
	\item{p1}{vector of probabilities corresponding to \code{x1}.
	}
  	\item{x2}{
		support vector of the second distribution, with strictly increasing elements.
	}
	\item{p2}{vector of probabilities corresponding to \code{x2}.
	}
}
\details{
The convolution is performed in C, looping through all paired sums, 
augmenting existing values or inserting them with an update of the 
corresponding probabilities.
}
\value{
A matrix with first column the new support vector and the second column the 
corresponding probability vector.
}



\examples{
x1 <- c(1,2,3.5)
p1 <- c(.2,.3,.5)
x2 <- c(0,2.3,3,4)
p2 <- c(.1,.3,.3,.3)


conv(x1,p1,x2,p2)
}


