\name{Opt}
\alias{Opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Opt(Style = c("European", "American", "Asian", "Binary", "AverageStrike", "Barrier", "Chooser", "Compound", "DeferredPayment", "ForeignEquity", "ForwardStart", "Gap", "HolderExtendible", "Ladder", "Lookback", "MOPM", "Perpetual", "Quotient", "Rainbow", "Shout", "SimpleChooser", "VarianceSwap"), Right = c("Call", "Put", "Other"), S0 = 50, ttm = 2, K = 52, Curr = "$", ContrSize = 100, SName = "A stock share", SSymbol = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Style}{
%%     ~~Describe \code{Style} here~~
}
  \item{Right}{
%%     ~~Describe \code{Right} here~~
}
  \item{S0}{
%%     ~~Describe \code{S0} here~~
}
  \item{ttm}{
%%     ~~Describe \code{ttm} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{Curr}{
%%     ~~Describe \code{Curr} here~~
}
  \item{ContrSize}{
%%     ~~Describe \code{ContrSize} here~~
}
  \item{SName}{
%%     ~~Describe \code{SName} here~~
}
  \item{SSymbol}{
%%     ~~Describe \code{SSymbol} here~~
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
function (Style = c("European", "American", "Asian", "Binary", 
    "AverageStrike", "Barrier", "Chooser", "Compound", "DeferredPayment", 
    "ForeignEquity", "ForwardStart", "Gap", "HolderExtendible", 
    "Ladder", "Lookback", "MOPM", "Perpetual", "Quotient", "Rainbow", 
    "Shout", "SimpleChooser", "VarianceSwap"), Right = c("Call", 
    "Put", "Other"), S0 = 50, ttm = 2, K = 52, Curr = "$", ContrSize = 100, 
    SName = "A stock share", SSymbol = "") 
{
    Style = match.arg(Style)
    s = list(Name = Style)
    s$Vanilla = (Style == "European" || Style == "American")
    s$Exotic = !s$Vanilla
    s$European = (Style == "European")
    s$American = (Style == "American")
    s$Asian = (Style == "Asian")
    s$Binary = (Style == "Binary")
    s$AverageStrike = (Style == "AverageStrike")
    s$Barrier = (Style == "Barrier")
    s$Chooser = (Style == "Chooser")
    s$Compound = (Style == "Compound")
    s$DP = (Style == "DeferredPayment")
    s$ForeignEquity = (Style == "ForeignEquity")
    s$ForwardStart = (Style == "ForwardStart")
    s$Gap = (Style == "Gap")
    s$HolderExtendible = (Style == "HolderExtendible")
    s$Ladder = (Style == "Ladder")
    s$Lookback = (Style == "Lookback")
    s$MOPM = (Style == "MOPM")
    s$Perpetual = (Style == "Perpetual")
    s$Quotient = (Style == "Quotient")
    s$Rainbow = (Style == "Rainbow")
    s$Shout = (Style == "Shout")
    s$SimpleChooser = (Style == "SimpleChooser")
    s$VarianceSwap = (Style == "VarianceSwap")
    Right = match.arg(Right)
    r = list(Name = Right)
    r$Call = (Right == "Call")
    r$Put = (Right == "Put")
    r$Other = (Right == "Other")
    r$SignCP = (r$Call * 2 - 1) * if (r$Other) 
        NA
    else 1
    o = list(S0 = S0, ttm = ttm, K = K, Style = s, Right = r, 
        Curr = Curr, ContrSize = ContrSize, SName = SName, SSymbol = SSymbol)
    class(o) = "Opt"
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
