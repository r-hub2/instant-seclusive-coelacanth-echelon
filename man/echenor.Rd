\name{echenor}
\alias{echenor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Echelon spatial scan statistic based on Normal model
}
\description{
The \code{echenor} function detects spatial clusters using the echelon spatial scan statistic with a Normal model.
}
\usage{
echenor(echelon.obj, val, weight = NULL, K = length(val)/2, Kmin = 2, n.sim = 99,
        cluster.type = "high", cluster.legend.pos = "bottomleft",
        dendrogram = TRUE, cluster.info = FALSE, coo = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{echelon.obj}{
An object of class \code{echelon}. For details, see \code{\link{echelon}}.
}
  \item{val}{
A numeric vector of observed values, which may be positive or negative. \code{NA} values are not allowed.
}
  \item{weight}{
A numeric vector of weighted values (must be positive). If \code{NULL} (the default), all weights are set to 1. \code{NA} values are not allowed.
}
  \item{K}{
Maximum cluster size. If \code{K} >= 2 (integer), the cluster size is limited to \code{K} regions.
}
  \item{Kmin}{
Minimum cluster size. Must be at least 2, due to the use of a permutation-based Monte Carlo test.
}
  \item{n.sim}{
The number of Monte Carlo replications used for significance testing of detected clusters. If set to 0, significance is not assessed.
}
  \item{cluster.type}{
A character string specifying the cluster type. If \code{"high"}, the detected clusters have high rates (hotspot). If \code{"low"}, the detected clusters have low rates (coldspot).
}
  \item{cluster.legend.pos}{
The location of the legend on the dendrogram. (See \code{\link{legend}} for details.)
}
  \item{dendrogram}{
Logical. If TRUE, draws an echelon dendrogram with the detected clusters.
}
  \item{cluster.info}{
Logical. If TRUE, returns detailed results of the detected clusters.
}
  \item{coo}{
An array of (x, y) coordinates for the region centroids to plot a cluster map.
}
  \item{\dots}{
Related to dendrogram drawing. (See the help for \code{\link{echelon}})
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{clusters}{Each detected cluster.}
\item{scanned.regions}{A region list of all scanning processes.}
\item{simulated.LLR}{Monte Carlo samples of the log-likelihood ratio.}
%% ...
}
\references{
[1] Kulldorff M, Huang L, and Konty K. (2009). A scan statistic for continuous data based on the normal probability model. \emph{International Journal of Health Geographics}, \strong{8}, 58.

[2] Huang L, Tiwari R, Zuo J, Kulldorff M, and Feuer E. (2009) Weighted normal spatial scan statistic for heterogeneous population data. \emph{Journal of the American Statistical Association}, \strong{104}, 886--898.
}
\author{
Fumio Ishioka
}
\note{
Typical values of \code{n.sim} are 99, 999, 9999, ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{echelon}} for the echelon analysis.

\code{\link{echepoi}} for cluster detection based on echelons using Poisson model.

\code{\link{echebin}} for cluster detection based on echelons using Binomial model.
}

\examples{
##Hotspot detection for predicting median housing values in Boston

#Load data and build a predictive model
library(spData)
data("boston")
gp0 <- lm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data = boston.c)
pred <- predict(gp0, newdata = boston.c, se.fit = TRUE)

#Extract predicted values and calculate reliability weights (inverse of standard error)
V <- gp0$fitted.values
W <- 1 / (pred$se.fit + 1e-6)

#Load geospatial information for Boston housing tracts
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg",
                           package="spData")[1])
boston.nb <- spdep::poly2nb(boston.tr)

#Hotspot detection based on Normal model
gp.echelon <- echelon(x = V, nb = boston.nb)
gp.clusters <- echenor(gp.echelon, val = V, weight = W, K = 40, cluster.info = TRUE, 
  main = "Hgih value clusters", ens = FALSE)

#Extract detected clusters
MLC <- gp.clusters$clusters[[1]]
Secondary1 <- gp.clusters$clusters[[2]]
Secondary2 <- gp.clusters$clusters[[3]]

#Assign colors to clusters for plotting
cluster.col <- rep(0, length(V))
cluster.col[MLC$regionsID] <- 2
cluster.col[Secondary1$regionsID] <- 3
cluster.col[Secondary2$regionsID] <- 4

#Plot detected high-value clusters on a simple map
plot(boston.tr$geom, col = cluster.col, 
  main = "Detected high value clusters")
legend("bottomleft",
  legend = c(
    paste("1- p-value:", MLC$p),
    paste("2- p-value:", Secondary1$p),
    paste("3- p-value:", Secondary2$p)
  ),
  text.col = c(2, 3, 4)
)

#Interactive map visualization using mapview
library(mapview)
boston.tr$cluster.col <- cluster.col
mapview(boston.tr, zcol = "cluster.col", 
  col.regions=c("white", "red", "orange", "yellow"), 
  label = "TOWN", legend=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ echelon analysis }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial scan statistic }% __ONLY ONE__ keyword per line
\keyword{ spatial cluster deteciotn }% __ONLY ONE__ keyword per line
