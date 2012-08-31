\name{ShorelineFireEMS}
\alias{ShorelineFireEMS}
\docType{data}
\title{
Shoreline Fire and EMS Turnout Times
}
\description{
This data set gives turnout response times for Fire and EMS (Emergency Medical Services) 
dispatch calls to the Shoreline, WA, Fire Department in 2006. The turnout time refers to 
time elapsed between the emergency call dispatch and the crew leaving the fire station,
or signaling that they are on their way while being on route already. The latter scenario
may explain the bimodal distribution character.
}
\usage{
data(ShorelineFireEMS)
}

\format{
  	A list of two sublists \code{$EMSTOT} and \code{$FireTOT}, 
	each with 4 vector components \code{$ST57}, \code{$ST63}, \code{$ST64}, and
	\code{$ST65} respectively, giving the turnout times (in seconds) (for EMS and Fire)
	at fire stations ST57, ST63, ST64, and ST65.
}
\note{These data sets are provided to illustrate usage of \code{ad.test} and \code{qn.test} 
and their combined versions in testing for performance equivalence across fire stations.
}
\source{
Thanks to Michael Henderson and the Fire Fighters
and Paramedics of the Shoreline Fire Department in Washington State.
}

\examples{
data(ShorelineFireEMS)
boxplot(ShorelineFireEMS$EMSTOT,xlab="Station", ylab="seconds",
	main="EMS Turnout Time")
boxplot(ShorelineFireEMS$FireTOT,xlab="Station", ylab="seconds",
	main="Fire Turnout Time")
}
\keyword{datasets}
