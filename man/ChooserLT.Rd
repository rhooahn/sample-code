\name{ChooserLT}
\alias{ChooserLT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChooserLT(o = OptPx(Opt(Style = "Chooser")), ttc = 0.5, IncBT = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{ttc}{
%%     ~~Describe \code{ttc} here~~
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
function (o = OptPx(Opt(Style = "Chooser")), ttc = 0.5, IncBT = FALSE) 
{
    stopifnot(is.OptPx(o), o$Style$Chooser, is.numeric(ttc))
    NSteps = o$NSteps
    p = o$p
    K = o$K
    step.ttc = ttc/o$ttm * NSteps
    S = with(o, S0 * d^(0:NSteps) * u^(NSteps:0))
    O.call = pmax((S - K), 0)
    O.put = pmax(-(S - K), 0)
    O.chooser = c()
    ReCalc.O_S.on.Prior.Time.Step = function(i) {
        if (i >= step.ttc + 1) {
            O.call <<- o$DF_dt * (p * O.call[-i - 1] + (1 - p) * 
                O.call[-1])
            O.put <<- o$DF_dt * (p * O.put[-i - 1] + (1 - p) * 
                O.put[-1])
            S <<- o$d * S[-i - 1]
            Payout.call = pmax((S - K), 0)
            Payout.put = pmax(-(S - K), 0)
            if (o$Style$American == TRUE) {
                O.call <<- pmax(O.call, Payout.call)
                O.put <<- pmax(O.put, Payout.put)
            }
            return(cbind(S, O.call, O.put))
        }
        else if (i == step.ttc) {
            S <<- o$d * S[-i - 1]
            O.chooser <<- o$DF_dt * (p * pmax(O.call[-i - 1], 
                O.put[-i - 1]) + (1 - p) * pmax(O.call[-1], O.put[-1]))
            Payout.call = pmax((S - K), 0)
            Payout.put = pmax(-(S - K), 0)
            if (o$Style$American == TRUE) {
                O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)
            }
            return(cbind(S, O.chooser))
        }
        else {
            S <<- o$d * S[-i - 1]
            O.chooser <<- o$DF_dt * (p * O.chooser[-i - 1] + 
                (1 - p) * O.chooser[-1])
            Payout.call = pmax((S - K), 0)
            Payout.put = pmax(-(S - K), 0)
            if (o$Style$American == TRUE) {
                O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)
            }
            return(cbind(S, O.chooser))
        }
    }
    BT = append(list(cbind(S, O.call, O.put)), sapply(NSteps:1, 
        ReCalc.O_S.on.Prior.Time.Step))
    o$ttc <- ttc
    o$IncBT <- IncBT
    o$PxLT = BT[[length(BT)]][[2]]
    if (IncBT) 
        o$BT = BT
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
