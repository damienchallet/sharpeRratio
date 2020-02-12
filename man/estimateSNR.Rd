\name{estimateSNR}
\alias{estimateSNR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A moment-free estimator of the Sharpe (signal-to-noise) ratio}

\description{
This function accepts a vector of price returns (or any possibly heavy-tailed data) and returns a list containing the moment-free estimate. 
}
\usage{
estimateSNR(x,numPerm=1000,nu=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{numPerm}{The basic assumption of the estimator is that the sample data are independent and indentically distributed. To improve the efficiency (precision) of the test, it is a good idea to average it over several random index permutations. You can use the default value (1000).}
  \item{nu}{If set to NA, the tail exponent of the data is obtained from fit to a Student t-distribution. If not NA, nu is used as the tail exponent of the Student t distribution in order to estimate the Sharpe ratio.}
}
\details{
The estimator first computes the cumulated sum of x (e.g. prices) and then counts the number of upper records of the price (starting from the first point);  an equivalent view point consists in computing the total drawdown duration.  For the sake of symmetry, the number of lower records (or total drawup duration) is also computed, and the estimator is simply the difference between the two.  Since it is an integer number, and provided that x is independently and identically distributed, the precision of the estimator may be improved by averaging the estimator on several random permutations of x, at the cost of computation time. This estimator assumes that the data is distributed according to Student t distribution.
}
\value{
\item{SNR}{The signal-to-noise ratio. To have something comparable with a t-statistics, multiply by sqrt(length(x)).}
\item{nu}{The fitted Student t-distribution tail exponent.}
\item{R0bar}{The number of upper records minus the number of lower records of the cumulated sum of x.}
\item{N}{The length of the vector x. It may be smaller than the input length if x contains NAs.}
}
\references{
Challet, D. (2017). Sharper asset ranking from total drawdown durations. Applied Mathematical Finance, 24(1), 1-22.
}
\author{
Damien Challet
}
\examples{
x <- rt(100,3)/sqrt(3)+0.05  #some Student-t distributed synthetic price log-returns
estimateSNR(x)    
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }% use one of  RShowDoc("KEYWORDS")
\keyword{ htest }% __ONLY ONE__ keyword per line
%@keywords ts htest
