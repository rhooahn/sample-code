\name{LookbackBS}
\alias{LookbackBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LookbackBS(o = OptPx(Opt(Style = "Lookback")), Smax = 50, Smin = 50, Type = c("Floating", "Fixed"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{Smax}{
%%     ~~Describe \code{Smax} here~~
}
  \item{Smin}{
%%     ~~Describe \code{Smin} here~~
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
function (o = OptPx(Opt(Style = "Lookback")), Smax = 50, Smin = 50, 
    Type = c("Floating", "Fixed")) 
{
    stopifnot(is.OptPx(o), o$Style$Name == "Lookback", is.numeric(Smax), 
        is.numeric(Smin), is.character(Type))
    Type = match.arg(Type)
    isFixed = switch(Type, Fixed = TRUE, Floating = FALSE)
    o$Smax = Smax
    o$Smin = Smin
    o$Type = Type
    fput.px <- function(S.max) {
        b1 = with(o, (log(S.max/S0) + (-r + q + vol^2/2) * ttm)/(vol * 
            sqrt(ttm)))
        b2 = with(o, b1 - vol * sqrt(ttm))
        b3 = with(o, (log(S.max/S0) + (r - q - vol^2/2) * ttm)/(vol * 
            sqrt(ttm)))
        Y2 = with(o, (2 * (r - q - vol^2/2) * log(S.max/S0))/vol^2)
        p_fl = with(o, S.max * exp(-r * ttm) * (pnorm(b1) - (vol^2 * 
            exp(Y2) * pnorm(-b3))/(2 * (r - q))) + (S0 * exp(-o$q * 
            ttm) * vol^2 * pnorm(-b2))/(2 * (r - q)) - S0 * exp(-q * 
            ttm) * pnorm(b2))
        return(p_fl)
    }
    fcall.px <- function(S.min) {
        a1 = with(o, (log(S0/S.min) + (r - q + vol^2/2) * ttm)/(vol * 
            sqrt(ttm)))
        a2 = with(o, a1 - vol * sqrt(ttm))
        a3 = with(o, (log(S0/S.min) + (-r + q + vol^2/2) * ttm)/(vol * 
            sqrt(ttm)))
        Y1 = with(o, (-2 * (r - q - vol^2/2) * log(S0/S.min))/vol^2)
        c_fl = with(o, S0 * exp(-q * ttm) * pnorm(a1) - (S0 * 
            exp(-q * ttm) * vol^2 * pnorm(-a1))/(2 * (r - q)) - 
            S.min * exp(-r * ttm) * (pnorm(a2) - (vol^2 * exp(Y1) * 
                pnorm(-a3))/(2 * (r - q))))
        return(c_fl)
    }
    if (o$Right$Name == "Put") {
        o$PxBS <- with(o, if (!isFixed) 
            fput.px(Smin)
        else fcall.px(min(Smin, K)) + K * exp(-r * ttm) - S0 * 
            exp(-q * ttm))
    }
    else if (o$Right$Name == "Call") {
        o$PxBS <- with(o, if (!isFixed) 
            fcall.px(Smax)
        else fput.px(max(Smax, K)) + S0 * exp(-q * ttm) - K * 
            exp(-r * ttm))
    }
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
