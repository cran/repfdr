\name{ztobins}
\alias{ztobins}
\title{
Binning of z-scores and estimation of the probabilities in each bin for the null and non-null states.
}
\description{
For each study, the function discretizes the z-scores into bins and estimates the probabilities in each bin for the null and non-null states using the function \code{\link{locfdr}}.
}
\usage{
ztobins(zmat, n.association.status = 3, n.bins = 120,nulltype = 0, ...)
}
\arguments{
  \item{zmat}{
Matrix of z-scores of the features (in rows) in each study (columns).
}
  \item{n.association.status}{
either 2 for no-association\\association or 3 for no-associtation\\negative-association\\positive-association.
}
  \item{n.bins}{
Number of bins in the discretization of the z-score axis.
}
  \item{nulltype}{
Type of null hypothesis assumed in estimating \code{f0(z)}. \code{0} (default) is the theoretical null N(0,1). This is an argument to pass to \code{locfdr}. For other values, see documentation of \code{\link{locfdr}}.
}
  \item{\dots}{
Arguments to pass to \code{locfdr}. See \code{\link{locfdr}} for details.
}
}
\details{
This function outputs the first two arguments to be input in the main function \code{\link{repfdr}}.
}
\value{
A list with:
  \item{pdf.binned.z}{
A 3-dimensional array which contains for each study (first dimension), the probabilities of a z-score to fall in the bin (second dimension), under each hypothesis status (third dimension). The third dimension can be of size 2 or 3, depending on the number of association states: if the association can be either null or only in one direction, the dimension is 2; if the association can be either null, or positive, or negative, the dimension is 3.  
}
  \item{binned.z.mat}{
A matrix of the bin numbers for each the z-scores (rows) in each study (columns).
}
}
\seealso{
\code{\link{repfdr}},\code{\link{locfdr}}
}
\examples{
data(zmat)

# three association states case (H in {-1,0,1}):
input.to.repfdr3 <- ztobins(zmat, 3, plot = TRUE, df = 15)
pbz    <- input.to.repfdr3$pdf.binned.z
bz     <- input.to.repfdr3$binned.z.mat

## Simulation: two association states case (H in {0,1}):

# data generation:
H <- hconfigs(n.studies= 3, n.association.status=2)
f <- c(0.895,0.005,0.005,0.02,0.005,0.02,0.02,0.03) 
cbind(H,f) # the simulation design
sum(f)     # all sum to 1?

m = 100000  # 100000 tests in each study
freq <- m*c(f,f,f)
Hvec <- rep(H,freq) 
set.seed(12)
simzmat  <- matrix(rnorm(n=3*m,mean=3*Hvec),nrow=m,ncol=3,byrow=FALSE)

# which of the tests are true replication\\association?
true.rep   <- rep((rowSums(abs(H)) > 1),m*f) 
true.assoc <- rep((rowSums(abs(H)) >= 1),m*f) 

input.to.repfdr <- ztobins(simzmat, 2, plot = TRUE, df = 15)
pbz    <- input.to.repfdr$pdf.binned.z
bz     <- input.to.repfdr$binned.z.mat
}
