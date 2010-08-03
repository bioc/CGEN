% --- Source file: snp.effects.Rd ---
\name{snp.effects}
\alias{snp.effects}
\title{Joint and Stratified Effects }
\description{
    Computes joint and stratified effects of the SNP and another variable based on a fitted model.}
\usage{ snp.effects(fit, var, var.levels=c(0, 1), method=NULL)}
\arguments{
   \item{fit}{ Return object from \code{\link{snp.logistic}} or \code{\link{snp.matched}}.
              If \code{fit} is the return object from \code{\link{snp.matched}}, then the \code{snp.vars}
              argument in \code{\link{snp.matched}} must consist of a single SNP.
              No default.}
   \item{var}{Name of the second variable to compute the effects for. This variable can be 
        a dummy variable, continuous variable, or a factor. Note that if this variable enters the model
      as both a main effect and interaction, then it must enter the model the same way as a main effect
     and interaction for the effects to be computed correctly.
     For example, if \code{var} is a factor as a main effect, then it also must 
    be a factor as an interaction. No default. }
   \item{var.levels}{(For continuous \code{var}) Vector of levels.
             First level is assumed to be the baseline level. 
            The default is c(0, 1).}
  \item{method}{Vector of values from "UML", "CML", "EB" or "CCL", "HCL", "CLR". The default is NULL.}
}
\details{The joint and stratified effects are computed for each method in \code{fit}.
         The stratified effects are the sub-group effect of the SNP stratified by
         \code{var} and the sub-group effect of \code{var} stratified by the SNP. \cr \cr
\bold{Definition of joint and stratified effects:} \cr
Consider the model: 
\deqn{logit(P(y=1)) = \alpha + \beta SNP + \gamma X + \delta SNP X.}{logit(P(y=1)) = alpha + beta*SNP + gamma*X + delta*SNP*X.}
Let 0 be the baseline for SNP and \eqn{x_0}{x_0} the baseline for X. Then the joint effect 
for SNP = s and X = x relative to SNP = 0 and X = \eqn{x_0}{x_0} is
\deqn{ \frac{\exp(\alpha + \beta s + \gamma x + \delta s x)}{\exp(\alpha + \gamma x_0)}}{exp(alpha + beta*s + gamma*x + delta*s*x)/exp(alpha + gamma*x_0)} 
The stratified effect of the SNP relative to SNP = 0 given X = x is
\deqn{ \frac{\exp(\alpha + \beta s + \gamma x + \delta s x)}{\exp(\alpha + \gamma x)}}{exp(alpha + beta*s + gamma*x + delta*s*x)/exp(alpha + gamma*x)} 
The stratified effect of \code{var} relative to X = x0 given SNP = s is
\deqn{ \frac{\exp(\alpha + \beta s + \gamma x + \delta s x)}{\exp(\alpha + \beta s)}}{exp(alpha + beta*s + gamma*x + delta*s*x)/exp(alpha + beta*s)} 
A convenient way to print the returned object to view the effects tables is with the function \code{\link{printEffects}}.
}

\value{
   If \code{fit} is of class \code{snp.logistic}, then the return object is a list of with names "UML", "CML", and "EB".
   If \code{fit} is of class \code{snp.matched}, then the return object is a list of with names "CLR", "CCL", and "HCL".
   Each sublist contains joint effects, stratified effects, standard errors and confidence intervals. 
   The sub-group effect of the SNP stratified by \code{var} is in the list "StratEffects", and the
   sub-group effect of \code{var} stratified by the SNP is in the list "StratEffects.2".
}
%\references{ }
%\author{ }
\seealso{ \code{\link{printEffects}} }
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 # Fit using a stratification variable
 fit <- snp.logistic(Xdata, "case.control", "BRCA.status",
                     main.vars=c("oral.years", "n.children"), 
                     int.vars=c("oral.years", "n.children"), 
                     strata.var="ethnic.group")

 # Compute the effects
 effects <- snp.effects(fit, "oral.years", var.levels=0:5) 
}
\keyword{ misc } 



