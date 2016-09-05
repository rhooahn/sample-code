#'@title       Chooser option valuation via Black-Scholes (BS) model 
#'@description Compute an exotic option that allow the holder decide the option 
#'              will be a call or put option at some predetermined future date.
#'              In a simple case, both put and call option are plain vanilla option.
#'      The value of the simple chooser option is \eqn{\max{C(S,K,t_1),P(S,K,t_2)}}.
#'              The plain vanilla option is calculated based on the BS model.
#'@author      Le You, Department of Statistics, Rice University, spring 2015
#'              
#'@param        o An object of class \code{OptPx}
#'@param        t1    The time to maturity of the call option, measured in years.
#'@param        t2    The time to maturity of the put option, measured in years.
#'@return       A list of class \code{SimpleChooserBS} consisting of the original \code{OptPx} object 
#'              and the option pricing parameters \code{t1}, \code{t2}, 
#'              as well as the computed price \code{PxBS}.
#'              
#'@references   
#'\itemize{
#' \item{Hull, John C.,\emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'      ISBN 978-0-13-345631-8. \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}}
#' \item{Huang Espen G., \emph{Option Pricing Formulas}, 2ed. 
#'     \url{http://down.cenet.org.cn/upfile/10/20083212958160.pdf} }
#' \item{Wee, Lim Tiong, MFE5010 \emph{Exotic Options,Notes for Lecture 4 Chooser option}.
#'     \url{http://www.stat.nus.edu.sg/~stalimtw/MFE5010/PDF/L4chooser.pdf} }
#' \item{Humphreys, Natalia A., ACTS 4302 Principles of Actuarial Models: Financial Economics.
#'     \emph{Lesson 14: All-or-nothing, Gap, Exchange and Chooser Options}. 
#'      \url{http://www.utdallas.edu/~natalia.humphreys/MFE AU14/MFE Lesson 14 slides.pdf}} }
#'              
#'@examples 
#'(o = ChooserBS())$PxBS
#'
#'o = Opt(Style='Chooser',Right='Other',S0=50, K=50)
#'(o = ChooserBS(OptPx(o, r=0.06, q=0.02, vol=0.2),9/12, 3/12))$PxBS
#'  
#'o = Opt(Style='Chooser',Right='Other',S0=50, K=50)
#'(o = ChooserBS (OptPx(o,r=0.08, q=0, vol=0.25),1/2, 1/4))$PxBS  
#'
#'o = Opt(Style='Chooser',Right='Other',S0=100, K=50)
#'(o = ChooserBS(OptPx(o,r=0.08, q=0.05, vol=0.3),1/2, 1/4))$PxBS
#'
#'@export 
#'
ChooserBS =function(o=OptPx(Opt(Style='Chooser')),t1=9/12,t2=3/12)
{ stopifnot(is.OptPx(o), o$Style$Chooser)
  d2 = (log(o$S0/o$K) + ((o$r-o$q) - o$vol ^ 2 / 2) * (t1)) / (o$vol * sqrt(t1))
  d1 = d2 + o$vol* sqrt(t1)
  d2n = ((log(o$S0/o$K) + (o$r-o$q) * t1 - o$vol^ 2 * t2 / 2) /(o$vol * sqrt(t2)))
  d1n = d2n + o$vol* sqrt(t2)
  o$PxBS =(o$S0 * exp (-o$q * t1) * pnorm(d1) - o$K * exp(-o$r * t1)* pnorm(d2)  + o$K * exp(-o$r * t1) * pnorm(-d2n) - o$S0 * exp (-o$q *t1) * pnorm(-d1n))
  return(o)
}



#' @title Chooser option valuation via Lattice Tree (LT) Model
#' @description Calculates the price of a Chooser option using a recombining binomial tree model. Has pricing capabilities for both simple European Chooser options as well as American Chooser Options, where exercise can occur any time as a call or put options.
#'
#' @author Richard Huang, Department of Statistics, Rice University, spring 2015
#' 
#' @param o The \code{OptPx} option object to price. 
#' @param ttc Time to choice: time (in years) until choice between call or put style option is required.
#' @param IncBT \code{TRUE/FALSE} Choice of including the lattice tree simulation in the output. 
#' Input \code{'FALSE'} yields faster computation and fewer calculated results to store in memory.
#' 
#' @details The \code{'American'} chooser option is interpreted as exercise of option being available at any point in time during the life of the option.
#' 
#' @return An original \code{OptPx} object with \code{PxLT} field as the price of the option and user-supplied \code{ttc}, 
#' \code{IncBT} parameters attached.
#'
#' @examples
#' (o = ChooserLT())$PxLT #Default Chooser option price. (See Ho pg 234 in references)
#' 
#' o = Opt('Eu', S0=100, ttm=1, K=100)
#' o = OptPx(o, r=0.10, q=0, vol=0.1, NSteps=500)
#' ChooserLT(o, ttc = .5, IncBT=T) 
#' 
#' #American Chooser, higher price than European equivalent
#' o = Opt('Am', S0=100, ttm=1, K=100)
#' o = OptPx(o, r=0.10, q=0, vol=0.1, NSteps=6)
#' ChooserLT(o,ttc=.5,IncBT=T) 
#' 
#' o = Opt('Eu', S0=50, ttm=1, K=50)
#' o = OptPx(o, r=0.05, q=0.02, vol=0.25, NSteps=1000)
#' ChooserLT(o, ttc = .75, IncBT=T) 
#' 
#' o = Opt('Eu', S0=50, ttm=1, K=50)
#' o = OptPx(o, r=0.05, q=0.5, vol=0.25, NSteps=100)
#' ChooserLT(o, ttc = .75, IncBT=F) 
#'  
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}, 
#' Thomas S.Y. Ho et al., \emph{The Oxford Guide to Financial Modeling : Applications for Capital Markets. . .}
#'
#' @export
#' 
ChooserLT = function(o=OptPx(Opt(Style='Chooser')), ttc = .5, IncBT=FALSE){ 
  stopifnot(is.OptPx(o), o$Style$Chooser, is.numeric(ttc)) #, (ttc/o$ttm * o$NSteps)%%1 == 0); 
  NSteps=o$NSteps; p=o$p; K=o$K
  
  step.ttc = ttc/o$ttm * NSteps
  S = with(o, S0*d^(0:NSteps)*u^(NSteps:0)) # vector of terminal stock prices, highest to lowest (@t=ttm)
  O.call = pmax((S - K), 0) # vector of terminal call option payouts (@t=ttm)
  O.put = pmax(-(S - K), 0) # vector of terminal put option payouts (@t=ttm)
  O.chooser = c()
  
  ReCalc.O_S.on.Prior.Time.Step = function(i) { 
    if (i >= step.ttc + 1) {
      O.call <<- o$DF_dt * (p*O.call[-i-1] + (1-p)*O.call[-1])  #prior call option prices (@time step=i-1)
      O.put <<- o$DF_dt * (p*O.put[-i-1] + (1-p)*O.put[-1])     #prior put option prices (@time step =i-1)
      S <<- o$d * S[-i-1]                   # prior stock prices (@time step=i-1)
      
      Payout.call = pmax((S - K), 0)   # call payout at time step i-1 (moving backward in time)
      Payout.put = pmax(-(S - K), 0)   # put payout at time step i-1 (moving backward in time)
      if (o$Style$American == TRUE) { 
        O.call <<- pmax(O.call, Payout.call)    #Test for American Chooser, max of payout values & expected values 
        O.put <<- pmax(O.put, Payout.put)
      }
      return(cbind(S, O.call, O.put))
    }
    else if (i == step.ttc) { # Time step where choice between put/call is made
      S <<- o$d * S[-i-1]                   # prior stock prices (@time step=i-1)
      O.chooser <<- o$DF_dt * (p*pmax(O.call[-i-1], O.put[-i-1]) + (1-p)*pmax(O.call[-1], O.put[-1]))  #prior chooser option prices (@time step=i-1)        
      Payout.call = pmax((S - K), 0)   
      Payout.put = pmax(-(S - K), 0)
      if (o$Style$American == TRUE) {
        O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)    #
      }
      return(cbind(S,O.chooser))
    }
    else { #Rolling back to current time
      S <<- o$d * S[-i-1]                   
      O.chooser <<- o$DF_dt * (p*O.chooser[-i-1] + (1-p)*O.chooser[-1])  
      Payout.call = pmax((S - K), 0)  
      Payout.put = pmax(-(S - K), 0)
      if (o$Style$American == TRUE) {
        O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)    #
      }
      return(cbind(S,O.chooser))
    }
  }
  
  BT = append(list(cbind(S, O.call, O.put)), sapply(NSteps:1, ReCalc.O_S.on.Prior.Time.Step)) #binomial tree
  
  #Return
  o$ttc <- ttc
  o$IncBT <- IncBT
  o$PxLT = BT[[length(BT)]][[2]]  # add ChoosterLT price
  if (IncBT) o$BT = BT
  return(o) 
}


