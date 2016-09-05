\name{BOPM_Eu}
\alias{BOPM_Eu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BOPM_Eu(o = OptPx())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
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
function (o = OptPx()) 
{
    stopifnot(is.OptPx(o), o$Style$European)
    with(o, {
        S = S0 * d^(NSteps:0) * u^(0:NSteps)
        O = pmax(o$Right$SignCP * (S - K), 0)
        csl = cumsum(log(c(1, 1:NSteps)))
        tmp = csl[NSteps + 1] - csl - csl[(NSteps + 1):1] + log(p) * 
            (0:NSteps) + log(1 - p) * (NSteps:0)
        o$PxBT = DF_ttm * sum(exp(tmp) * O)
        return(o)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
