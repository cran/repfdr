\name{piem}
\alias{piem}
\alias{Pi}
\title{
Estimation of the prior probabilities for each association status vector.
}
\description{
The function calls an expectation-maximization (EM) algorithm to estimate the prior probabilities of each association status vector. It is also used internally in \code{\link{repfdr}}.
}
\usage{
piem(pdf.binned.z, binned.z.mat, control = em.control())}
\arguments{
  \item{pdf.binned.z}{Same input as in \code{\link{repfdr}}.
 A 3-dimensional array which contains for each study (first dimension), the probabilities of a z-score to fall in the bin (second dimension), under each hypothesis status (third dimension). The third dimension can be of size 2 or 3, depending on the number of association states: if the association can be either null or only in one direction, the dimension is 2; if the association can be either null, or positive, or negative, the dimension is 3.  
Element \code{[[1]]} in the output of \code{\link{ztobins}}. 
}
  \item{binned.z.mat}{
  Same input as in \code{\link{repfdr}}.
A matrix of the bin numbers for each the z-scores (rows) in each study (columns).
Element \code{[[2]]} in the output of \code{\link{ztobins}}.
}
  \item{control}{
List of control parameters to pass to the EM algorithm. See \code{\link{em.control}}.
}
}
\details{
The implementation of the EM algorithm is in C, and allows paralel processing. By default, the software automatically detects the number of available processing threads. See \code{\link{em.control}} for the option of providing the number of threads to use, as well as for the additional control parameters. 
}
\value{
  \item{all.iterations}{Matrix with number of columns equal to  the number of EM iterations, and each column is the estimated probability distribution of the vector of association status.}
  \item{last.iteration}{Matrix of the vectors of association status along with the column vector of the last EM iteration, which contains the estimated probabilities of the vectors of association status.}
}
\references{
Heller, Ruth, and Daniel Yekutieli. "Replicability analysis for Genome-wide Association studies." \emph{arXiv preprint arXiv:1209.2829} (2012).
}
\author{
\code{C} implementation by Shachar Kaufman.
}
\seealso{
\code{\link{repfdr}}
}
\examples{
\dontrun{
data(zmat)
input.to.piem <- ztobins(zmat, 3, plot = TRUE)
pbz    <- input.to.piem$pdf.binned.z
bz     <- input.to.piem$binned.z.mat

output <- piem(pbz, bz, control = em.control(verbose = TRUE,nr.threads = 2))
# converge after 209 iterations.
# iterations are printed
# Set nr.threads = 0 to optimal use of threads.

Pi <- output$last.iteration
}
}
