% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{moc.gabk}
\alias{moc.gabk}
\title{Perform the Multi-Objective Clustering Algorithm Guided by a-Priori Biological Knowledge (MOC-GaPBK)}
\usage{
moc.gabk(
  dmatrix1,
  dmatrix2,
  num_k,
  generation = 50,
  pop_size = 10,
  rat_cross = 0.8,
  rat_muta = 0.01,
  tour_size = 2,
  neighborhood = 0.1,
  local_search = FALSE,
  cores = 2
)
}
\arguments{
\item{dmatrix1}{A distance matrix. It should have the same dimensions that dmatrix2. It is mandatory.}

\item{dmatrix2}{A distance matrix. It should have the same dimensions that dmatrix1. It is mandatory.}

\item{num_k}{The number k of groups represented by medoids in each individual. It is mandatory.}

\item{generation}{Number of generations to be performed by MOC-GaPBK. By default 50.}

\item{pop_size}{Size of population. By default 10.}

\item{rat_cross}{Probability of crossover. By default 0.80.}

\item{rat_muta}{Probability of mutation. By default 0.01.}

\item{tour_size}{Size of tournament. By default 2.}

\item{neighborhood}{Percentage of neighborhood. A real value between 0 and 1. It is computed as neighborhood*pop_size to determine the size of neighborhood. By default 0.10.}

\item{local_search}{A boolean value indicating whether the local searches procedures (PR and PLS) are computed. By default \emph{FALSE}.}

\item{cores}{Number of cores to be used to compute the local searches procedures. By default 2.}
}
\value{
\item{population}{The population of medoids including the objective functions values and order by Pareto ranking and crowding distance values.}
\item{matrix.solutions}{A matrix with results of clustering. Each column represents a clustering solution available in Pareto front.}
\item{clustering}{A list containing named vectors of integers from 1:k representing the cluster to which each object is assigned.}
}
\description{
This function receives two distance matrices and it performs the MOC-GaPBK.
}
\details{
MOC-GaPBK is a method proposes by Parraga-Alava, J. et. al. 2018. It carries out the discovery of clusters using NSGA-II algorithm along with Path-Relinking (PR) and Pareto Local Search (PLS) as intensification and diversification strategies, respectively. The algorithm uses as objective functions two versions of the Xie-Beni validity index, i.e., a version for each distance matrix (dmatrix1, dmatrix2). More details about this compute can be found in: <https://doi.org/10.1186/s13040-018-0178-4>. MOC-GaPBK yield a set of the best clustering solutions from a multi-objective point of views.
}
\examples{

##Generates a data matrix of dimension 50X20

library("amap")
library("moc.gapbk")

x <- matrix(runif(50 * 20, min = -5, max = 10), nrow = 50, ncol = 20)

##Compute two distance matrices

dmatrix1<- as.matrix(amap::Dist(x, method = "euclidean"))
dmatrix2<- as.matrix(amap::Dist(x, method = "correlation"))

##Performs MOC-GaPBK with 5 cluster

example<-moc.gabk(dmatrix1, dmatrix2, 5)

example$population
example$matrix.solutions
example$clustering

}
\references{
J. Parraga-Alava, M. Dorn, M. Inostroza-Ponta (2018). \emph{A multi-objective gene clustering algorithm guided by apriori biological knowledge with intensification and diversification strategies}. BioData Mining. 11(1) 1-16.

K. Deb, A. Pratap, S. Agarwal, T. Meyarivan (2002). \emph{A fast and elitist multiobjective genetic algorithm: NSGA-II}. IEEE Transactions on Evolutionary Computation, 6(2) 182-197.

F. Glover (1997). \emph{Tabu Search and Adaptive Memory Programming - Advances, Applications and Challenges}. Interfaces in Computer Science and Operations Research: Advances in Metaheuristics, Optimization, and Stochastic Modeling Technologies. 1-75.

J. Dubois-Lacoste, M. Lopez-Ibanez, Stutzle, T. (2015). \emph{Anytime Pareto local search}. European Journal of Operational Research, 243(2) 369-385.
}
\author{
Jorge Parraga-Alava, Marcio Dorn, Mario Inostroza-Ponta
}
\keyword{Crowding}
\keyword{GaPBK}
\keyword{NSGA-II}
\keyword{Pareto}
\keyword{Xie-Beni}
\keyword{distance}
\keyword{index}
\keyword{ranking}
\keyword{validity}
