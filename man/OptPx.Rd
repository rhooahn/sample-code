\name{OptPx}
\alias{OptPx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OptPx(o = Opt(), r = 0.05, q = 0, rf = 0, vol = 0.3, NSteps = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{rf}{
%%     ~~Describe \code{rf} here~~
}
  \item{vol}{
%%     ~~Describe \code{vol} here~~
}
  \item{NSteps}{
%%     ~~Describe \code{NSteps} here~~
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
function (o = Opt(), r = 0.05, q = 0, rf = 0, vol = 0.3, NSteps = 2) 
{
    stopifnot(is.Opt(o), is.numeric(r), r > 0, r < 1, is.numeric(vol), 
        is.numeric(q), is.numeric(rf))
    dt = o$ttm/NSteps
    u = exp(vol * sqrt(dt))
    d = 1/u
    SYld = r - q - rf
    a = exp(SYld * dt)
    o$r = r
    o$q = q
    o$rf = rf
    o$vol = vol
    o$NSteps = NSteps
    o$u = u
    o$d = d
    o$dt = dt
    o$a = a
    o$p = p = (a - d)/(u - d)
    o$SYld = SYld
    o$DF_ttm = exp(-r * o$ttm)
    o$DF_dt = exp(-r * dt)
    class(o) = c(class(o), "OptPx")
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
