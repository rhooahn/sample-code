\name{ChooserBS}
\alias{ChooserBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChooserBS(o = OptPx(Opt(Style = "Chooser")), t1 = 9/12, t2 = 3/12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{t1}{
%%     ~~Describe \code{t1} here~~
}
  \item{t2}{
%%     ~~Describe \code{t2} here~~
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
function (o = OptPx(Opt(Style = "Chooser")), t1 = 9/12, t2 = 3/12) 
{
    stopifnot(is.OptPx(o), o$Style$Chooser)
    d2 = (log(o$S0/o$K) + ((o$r - o$q) - o$vol^2/2) * (t1))/(o$vol * 
        sqrt(t1))
    d1 = d2 + o$vol * sqrt(t1)
    d2n = ((log(o$S0/o$K) + (o$r - o$q) * t1 - o$vol^2 * t2/2)/(o$vol * 
        sqrt(t2)))
    d1n = d2n + o$vol * sqrt(t2)
    o$PxBS = (o$S0 * exp(-o$q * t1) * pnorm(d1) - o$K * exp(-o$r * 
        t1) * pnorm(d2) + o$K * exp(-o$r * t1) * pnorm(-d2n) - 
        o$S0 * exp(-o$q * t1) * pnorm(-d1n))
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
