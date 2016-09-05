\name{BOPM}
\alias{BOPM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BOPM(o = OptPx(), IncBT = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{IncBT}{
%%     ~~Describe \code{IncBT} here~~
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
function (o = OptPx(), IncBT = TRUE) 
{
    stopifnot(is.OptPx(o), o$Style$Vanilla)
    NSteps = o$NSteps
    p = o$p
    K = o$K
    if (o$Style$European && !IncBT) 
        return(BOPM_Eu(o))
    else {
        S = with(o, S0 * d^(0:NSteps) * u^(NSteps:0))
        O = pmax(o$Right$SignCP * (S - K), 0)
        RecalcOSonPriorTimeStep = function(i) {
            O <<- o$DF_dt * (p * O[-i - 1] + (1 - p) * O[-1])
            S <<- o$d * S[-i - 1]
            Payout = pmax(o$Right$SignCP * (S - K), 0)
            if (o$Style$American) 
                O <<- pmax(O, Payout)
            return(cbind(S, O))
        }
        BT = append(list(cbind(S, O)), sapply(NSteps:1, RecalcOSonPriorTimeStep))
        o$PxBT = BT[[length(BT)]][[2]]
        if (IncBT) 
            o$BT = BT
        return(o)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
