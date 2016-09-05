\name{LookbackMC}
\alias{LookbackMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LookbackMC(o = OptPx(Opt(Style = "Lookback"), r = 0.05, q = 0, vol = 0.3, ), NPaths = 10000, div = 1000, Type = c("Floating", "Fixed"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{NPaths}{
%%     ~~Describe \code{NPaths} here~~
}
  \item{div}{
%%     ~~Describe \code{div} here~~
}
  \item{Type}{
%%     ~~Describe \code{Type} here~~
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
function (o = OptPx(Opt(Style = "Lookback"), r = 0.05, q = 0, 
    vol = 0.3, ), NPaths = 10000, div = 1000, Type = c("Floating", 
    "Fixed")) 
{
    stopifnot(is(o, "OptPx"), o$Style$Name == "Lookback", is.character(Type), 
        is.numeric(NPaths), is.numeric(div))
    o.class = class(o)
    Type = match.arg(Type)
    isFixed = switch(Type, Fixed = T, Floating = F)
    isFloating = !isFixed
    dt = o$ttm/div
    payoffs = replicate(NPaths, {
        ds.s = (o$r - o$q) * dt + o$vol * sqrt(dt) * rnorm(div)
        st = c(o$S0, o$S0 * cumprod(1 + ds.s))
        smax = max(st)
        smin = min(st)
        if (isFloating) {
            K = smax * o$Right$Put + smin * o$Right$Call
            ST = st[length(st)]
        }
        else if (isFixed) {
            K = o$K
            ST = smin * o$Right$Put + smax * o$Right$Call
        }
        payoff = exp(-o$r * o$ttm) * o$Right$SignCP * (max(ST - 
            K))
    })
    p = mean(payoffs)
    o$isFixed = isFixed
    o$isFloating = !isFixed
    o$NPaths = NPaths
    o$PxMC = p
    class(o) = o.class
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
