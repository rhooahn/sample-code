\name{BS_Simple}
\alias{BS_Simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BS_Simple(S0 = 42, K = 40, r = 0.1, q = 0, ttm = 0.5, vol = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S0}{
%%     ~~Describe \code{S0} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{ttm}{
%%     ~~Describe \code{ttm} here~~
}
  \item{vol}{
%%     ~~Describe \code{vol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (S0 = 42, K = 40, r = 0.1, q = 0, ttm = 0.5, vol = 0.2) 
{
    stopifnot(is.numeric(S0), S0 >= 0, is.numeric(K), K > 0, 
        is.numeric(r), r >= 0 && r <= 1, is.numeric(q), q >= 
            0 && q <= 1, is.numeric(ttm), ttm > 0, is.numeric(vol), 
        vol > 0)
    d1 = (log(S0/K) + (r - q + (vol^2)/2) * ttm)/(vol * sqrt(ttm))
    d2 = d1 - vol * sqrt(ttm)
    Nd1 = pnorm(d1)
    Nd2 = pnorm(d2)
    Call = S0 * exp(-q * ttm) * Nd1 - K * exp(-r * ttm) * Nd2
    Put = K * exp(-r * ttm) * pnorm(-d2) - S0 * exp(-q * ttm) * 
        pnorm(-d1)
    Px = list(Call = Call, Put = Put)
    BS = list(d1 = d1, d2 = d2, Nd1 = Nd1, Nd2 = Nd2, Px = Px)
    return(BS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
