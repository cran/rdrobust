\name{rdbwselect}
\alias{rdbwselect}
\alias{print.rdbwselect}
\alias{summary.rdbwselect}
\alias{print.summary.rdbwselect}
\title{Bandwidth selection procedures for local-polynomial RD estimators}
\description{
    \code{\link{rdbwselect}} implements bandwidth selection procedures for local-polynomial RD estimators proposed in Calonico, Cattaneo and Titiunik (2014a). It also computes bandwith selection procedures proposed in Ludwig and Miller (2007) and Imbens and Kalyanaraman (2012).

    Calonico, Cattaneo and Titiunik (2014b) provides an introduction to this command. Additional details for conventional approaches to conduct inference in the RD design can be found in Imbens and Lemieux (2008), Lee and Lemieux (2010), Dinardo and Lee (2011), and references therein.

    For data-driven RD plots use the command \code{\link{rdplot}}, which implements the results in Calonico, Cattaneo and Titiunik (2014d).

    A companion \code{Stata} package is described in Calonico, Cattaneo and Titiunik (2014c).
}
\usage{rdbwselect(y, x, data, subset = NULL, c = 0, p = 1, q = 2, deriv = 0,
              rho = NULL, kernel = "tri", bwselect = "CCT", scaleregul = 1,
              delta = 0.5, cvgrid_min = NULL, cvgrid_max = NULL, 
              cvgrid_length = NULL, cvplot = FALSE, vce = "nn", matches = 3, 
              all = FALSE, precalc = TRUE, model = FALSE, frame = FALSE)}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the running variable (a.k.a. score or forcing variable).}  
  \item{data}{an optional data frame.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{c}{specifies the RD cutoff in \code{x}; default is \code{c = 0}.}
    \item{p}{specifies the order of the local-polynomial used to construct the point-estimator; default is \code{p = 1} (local linear regression).}
  \item{q}{specifies the order of the local-polynomial used to construct the bias-correction; default is \code{q = 2} (local quadratic regression).}
  \item{deriv}{specifies the order of the derivative of the regression function to be estimated; default is \code{deriv = 0} (Sharp RD, or Fuzzy RD if \code{fuzzy} is also specified). Setting it equal to 1 results in estimation of a Kink RD design (or Fuzzy Kink RD if \code{fuzzy} is also specified).}  
    \item{rho}{if specified, sets the pilot bandwidth \code{b} equal to \code{h/rho}, where \code{h} is computed using the method and options chosen below.}
  \item{kernel}{is the kernel function used to construct the local-polynomial estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
  \item{bwselect}{selects the bandwidth selection procedure to be used. By default it computes both \code{h} and \code{b}, unless \code{rho} is specified, in which case it only computes \code{h} and sets \code{b} = \code{h/rho}. Options are:
  
  \code{CCT} for bandwidth selector proposed by Calonico, Cattaneo and Titiunik (2014a)  (default option).
  
  \code{IK} for bandwidth selector proposed by Imbens and Kalyanaraman (2012) (only available for Sharp RD design).
  
  \code{CV} for cross-validation method proposded by Ludwig and Miller (2007) (only available for Sharp RD design).}
  \item{scaleregul}{specifies scaling factor for the regularization terms of \code{CCT} and \code{IK} bandwidth selectors. Setting \code{scaleregul = 0} removes the regularization term from the bandwidth selectors; default is \code{scaleregul = 1}.}
    \item{delta}{sets the quantile that defines the sample used in the cross-validation procedure. This option is used only if \code{bwselect(CV)} is specified; default is \code{delta = 0.5}, that is, the median of the control and treated samples.}
  \item{cvgrid_min}{sets the minimum value of the bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified.}
  \item{cvgrid_max}{sets the maximum value of the bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified.}
  \item{cvgrid_length}{sets the bin length of the (evenly-spaced) bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified.}
  \item{cvplot}{generates a graph of the CV objective function. This option is used only if \code{bwselect = "CV"} is specified.}
  \item{vce}{specifies the procedure used to compute the variance-covariance matrix estimator. This option is used only if \code{CCT} or \code{IK} bankdwith procedures are employed. Options are:

  \code{nn} for nearest-neighbor matches residuals using \code{matches} number of matches. This is the default option (with \code{matches = 3}, see below).

  \code{resid} for estimated plug-in residuals using \code{h} bandwidth.}
  \item{matches}{specifies the number of matches in the nearest-neighbor based variance-covariance matrix estimator. This options is used only when nearest-neighbor matches residuals are employed; default is \code{matches = 3}.}
  
\item{all}{if specified, \code{\link{rdbwselect}} reports  three different procedures: 

  \code{CCT} for bandwidth selector proposed by Calonico, Cattaneo and Titiunik (2014a).
  
  \code{IK} for bandwidth selector proposed by Imbens and Kalyanaraman (2012).
  
  \code{CV} for cross-validation method proposed by Ludwig and Miller (2007).}
  
  \item{model}{logical. If \code{TRUE}, the model object will be returned.}
  \item{frame}{logical. If \code{TRUE}, the data frame used  in model fitting will be returned.}

  \item{precalc}{internal option.}
}
\value{
  \item{bws}{matrix containing the estimated bandwidths for each selected procedure.}
  \item{bwselect}{bandwidth selection procedure employed.}
  \item{kernel}{kernel function used to construct the local-polynomial estimator(s).}
  \item{p}{order of the local-polynomial used to construct the point-estimator.}
  \item{q}{order of the local-polynomial used to construct the bias-correction estimator.}
}
\references{
        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014a. "Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs." \emph{Econometrica}, forthcoming.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.
        
        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014b. "rdrobust: An R Package for Robust Inference in Regression-Discontinuity Designs." University of Michigan, Department of Economics. \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Rpkg.pdf}.

    Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014c. "Robust Data-Driven Inference in the Regression-Discontinuity Design." University of Michigan, Department of Economics.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}.
    
     Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014d. "Optimal Data-Driven Regression Discontinuity Plots." University of Michigan, Department of Economics.   \url{http://www-personal.umich.edu/~cattaneo/papers/RD-rdplot.pdf}.

    Dinardo, J., and D. S. Lee. 2011. "Program Evaluation and Research Designs." In \emph{Handbook of Labor Economics}, ed. O. Ashenfelter and D. Card, vol. 4A, 463-536. Elsevier Science B.V.

    Imbens, G., and T. Lemieux. 2008. "Regression Discontinuity Designs: A Guide to Practice." \emph{Journal of Econometrics} 142(2): 615-635.

    Imbens, G. W., and K. Kalyanaraman. 2012. "Optimal Bandwidth Choice for the Regression Discontinuity Estimator." \emph{Review of Economic Studies} 79(3): 933-959.

    Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in Economics." \emph{Journal of Economic Literature} 48(2): 281-355.

    Ludwig, J., and D. L. Miller. 2007. "Does Head Start Improve Children's Life Chances? Evidence from a Regression Discontinuity Design." \emph{Quarterly Journal of Economics} 122(1): 159-208.
    
}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Rocio Titiunik, University of Michigan, Ann Arbor, MI.  \email{titiunik@umich.edu}.
}
\examples{
x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdbwselect(y,x)
}
\seealso{
  \code{\link{rdrobust}}, \code{\link{rdplot}}
}