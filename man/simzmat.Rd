\name{simzmat}
\alias{simzmat}
\alias{Hmat}
\docType{data}
\title{
Simulated data set 
}
\description{
A simulated data set from three studie, with 10000 "features" in each study, each of which yielded a z-score. The data comprises 10000x3 z-scores, along with a matrix of size 10000x3 of indicators of whether each z-score belongs to a non-null hypothesis for the feature in the study (1) or to a null hypothesis for the feature in the study (0).  
}
\usage{data(simzmat)}
\format{
\code{simzmat} is a matrix of 10000 z-scores (in rows) in each of the 3 studies (columns). 
\code{Hmat} is a matrix of 10000 rows, each row a vector of the true association status from which the z-scores in the same row in \code{simzmat} was generated. Specifically, for a zero entry in Hmat the corresponding z-score in \code{simzmat} was generated from the standard normal distribution, and for a unit entry in Hmat the corresponding z-score in \code{simzmat} was generated from the normal distribution with mean 3 and variance one.
}
\examples{
data(simzmat)

# the generation code:
\dontrun{
require(repfdr)

H <- hconfigs(n.studies= 3, n.association.status=2)
# frequencies for the association status vectors:
f <- c(0.895,0.005,0.005,0.02,0.005,0.02,0.02,0.03) 
cbind(H,f)  # the simulation design
sum(f)      # all sum to 1?
m = 10000  # 10000 tests in each study

# true states:
Hmat <- matrix(rep(x = H, times = m*cbind(f,f,f)),ncol=3) 

# simulated z-scores:
set.seed(12)
simzmat  <- matrix(rnorm(n=3*m,mean=Hmat*3),nrow=m,ncol=3)
rm(f,m,H)

}
}
\keyword{datasets}
\keyword{Replicability}
\keyword{GWAS}

