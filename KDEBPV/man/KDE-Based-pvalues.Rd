\name{KDEPV}
\alias{KDEPV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel Density Estimator (KDE) Based P-values
}
\description{
Computes p-values corresponding to KDE based Studentized t-statistics.
}
\usage{
KDEPV(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x: A data matrix
}
\item{y}{
y : A data matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns Benjamini and Hochberg adjusted p-values and the corresponding respective hypotheses in a matrix form.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Santu Ghosh (sghosh@augusta.edu)
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
  %}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
