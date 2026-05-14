#' Multi-Objective Clustering Guided by a-Priori Biological Knowledge (MOC-GaPBK)
#'
#' Performs the MOC-GaPBK algorithm proposed by Parraga-Alava and others
#' (2018). It receives two distance matrices and returns a set of
#' non-dominated clustering solutions.
#'
#' @param dmatrix1 A square distance matrix. Must have the same
#'   dimensions as \code{dmatrix2}.
#' @param dmatrix2 A square distance matrix. Must have the same
#'   dimensions as \code{dmatrix1}. Typically encodes a-priori
#'   biological knowledge.
#' @param num_k The number \eqn{k} of clusters represented by medoids in
#'   each individual. Must be greater than 1.
#' @param generation Number of generations to be performed. Default 50.
#' @param pop_size Size of the population. Default 10.
#' @param rat_cross Probability of crossover. Default 0.80.
#' @param rat_muta Probability of mutation. Default 0.01.
#' @param tour_size Size of the tournament for parent selection. Default 2.
#' @param neighborhood Percentage of neighborhood used by Pareto Local
#'   Search. A real value between 0 and 1. The neighborhood size is
#'   computed as \code{neighborhood * num_objects}. Default 0.10.
#' @param local_search Logical. If \code{TRUE}, Path-Relinking (PR) and
#'   Pareto Local Search (PLS) are applied as intensification and
#'   diversification strategies. Default \code{FALSE}.
#' @param cores Number of cores used by Path-Relinking. Default 2.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{population}}{A data frame containing the final
#'     population of medoids together with the values of the two
#'     objective functions, the Pareto ranking and the crowding
#'     distance, ordered accordingly.}
#'   \item{\code{matrix.solutions}}{A data frame whose columns are
#'     clustering solutions on the Pareto front. Each row corresponds
#'     to an object and each cell to its assigned cluster.}
#'   \item{\code{clustering}}{A list of named integer vectors. Element
#'     \code{i} is the partition produced by the \eqn{i}-th solution
#'     on the Pareto front.}
#' }
#'
#' @details
#' MOC-GaPBK couples NSGA-II with Path-Relinking and Pareto Local Search.
#' Two versions of the Xie-Beni validity index are used as objectives,
#' one per distance matrix.
#'
#' @author Jorge Parraga-Alava, Marcio Dorn, Mario Inostroza-Ponta
#'
#' @references
#' J. Parraga-Alava, M. Dorn, M. Inostroza-Ponta (2018). \emph{A
#' multi-objective gene clustering algorithm guided by apriori
#' biological knowledge with intensification and diversification
#' strategies}. BioData Mining. 11(1) 1-16.
#' \doi{10.1186/s13040-018-0178-4}.
#'
#' K. Deb, A. Pratap, S. Agarwal, T. Meyarivan (2002). \emph{A fast and
#' elitist multiobjective genetic algorithm: NSGA-II}. IEEE Transactions
#' on Evolutionary Computation, 6(2) 182-197.
#'
#' F. Glover (1997). \emph{Tabu Search and Adaptive Memory Programming -
#' Advances, Applications and Challenges}. Interfaces in Computer
#' Science and Operations Research. 1-75.
#'
#' J. Dubois-Lacoste, M. Lopez-Ibanez, T. Stutzle (2015). \emph{Anytime
#' Pareto local search}. European Journal of Operational Research,
#' 243(2) 369-385.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(stats::runif(50 * 20, min = -5, max = 10),
#'             nrow = 50, ncol = 20)
#'
#' # Two distance matrices from base R; in real applications dmatrix2
#' # typically encodes a-priori biological knowledge (e.g. GO semantic
#' # similarity). See vignette("moc-gapbk-intro") for examples using
#' # amap::Dist() with correlation-based distances.
#' dmatrix1 <- as.matrix(stats::dist(x, method = "euclidean"))
#' dmatrix2 <- as.matrix(stats::dist(x, method = "manhattan"))
#'
#' res <- moc.gapbk(dmatrix1, dmatrix2, num_k = 3,
#'                  generation = 5, pop_size = 6)
#'
#' head(res$matrix.solutions)
#'
#' @export
moc.gapbk <- function(dmatrix1, dmatrix2, num_k,
                      generation = 50, pop_size = 10,
                      rat_cross = 0.80, rat_muta = 0.01,
                      tour_size = 2, neighborhood = 0.10,
                      local_search = FALSE, cores = 2) {

  if (!(nrow(dmatrix1) == ncol(dmatrix1) &&
        nrow(dmatrix2) == ncol(dmatrix2))) {
    stop("'dmatrix1' or 'dmatrix2' must have the same number of rows and columns")
  }
  if (!(nrow(dmatrix1) == nrow(dmatrix2) &&
        ncol(dmatrix1) == ncol(dmatrix2))) {
    stop("Both matrices must have equal number of rows and columns")
  }
  if (missing(num_k)) stop("'num_k' must be indicated")
  if (num_k <= 1) stop("'num_k' should be more than 1")
  if (neighborhood < 0 || neighborhood > 1) {
    stop("'neighborhood' should be between 0 and 1")
  }

  num_objects <- nrow(dmatrix1)
  number.objectives <- 2

  # Performance: precompute squared distance matrices once. This is the
  # only quantity that calculate.objective.functions() needs, and the
  # original implementation was squaring elements element-by-element in
  # an inner for-loop. Caching here turns O(generation * pop_size * n)
  # squarings into a single O(n^2) prelude.
  attr(dmatrix1, "dist_sq") <- dmatrix1 * dmatrix1
  attr(dmatrix2, "dist_sq") <- dmatrix2 * dmatrix2

  population.P <- generate.initial.population(nrow(dmatrix1), num_k, pop_size)
  g <- 1

  while (g <= generation) {
    # Population P
    verify.singletons.P <- singletons.repair(population.P, dmatrix1, dmatrix2,
                                             num_objects)
    table.groups.P <- verify.singletons.P$groups
    population.P <- as.matrix(verify.singletons.P$population)
    population.P <- calculate.ranking.crowding(pop_size, population.P,
                                               table.groups.P,
                                               number.objectives, FALSE,
                                               dmatrix1, dmatrix2, num_k)

    mating.pool <- nsga2R::tournamentSelection(population.P, pop_size, tour_size)
    population.Q <- t(sapply(seq_len(pop_size),
                             function(u) array(rep(0, num_k))))

    crossover <- generate.crossover.k.points(pop_size, num_k, mating.pool,
                                             population.Q, rat_cross)
    population.Q <- verify.feasibility(crossover, population.Q, num_objects)

    mutation <- generate.mutation.random.controller(pop_size, num_k,
                                                    population.Q, rat_muta,
                                                    num_objects)
    population.Q <- verify.feasibility(mutation, population.Q, num_objects)

    # Population Q
    verify.singletons.Q <- singletons.repair(population.Q, dmatrix1, dmatrix2,
                                             num_objects)
    table.groups.Q <- verify.singletons.Q$groups
    population.Q <- as.matrix(verify.singletons.Q$population)
    population.Q <- calculate.ranking.crowding(pop_size, population.Q,
                                               table.groups.Q,
                                               number.objectives, FALSE,
                                               dmatrix1, dmatrix2, num_k)

    # Population R
    population.R <- rbind(population.P, population.Q)
    rownames(population.R) <- seq_len(nrow(population.R))
    population.R <- population.R[, -((num_k + 1):(num_k + number.objectives + 2))]

    table.groups.R <- generate.groups(nrow(population.R), population.R,
                                      dmatrix1, dmatrix2)
    population.R <- calculate.ranking.crowding(pop_size * 2, population.R,
                                               table.groups.R,
                                               number.objectives, FALSE,
                                               dmatrix1, dmatrix2, num_k)
    population.R <- as.data.frame(population.R)

    population.pareto <- subset(population.R, population.R$paretoranking == "1")
    population.pareto <- population.pareto[
      !duplicated(population.pareto[, seq_len(num_k)]),
    ]
    population.pareto <- population.R[
      !duplicated(population.R[, seq_len(num_k)]),
    ]
    table.groups.Pareto <- generate.groups(
      nrow(population.pareto),
      as.matrix(population.pareto[, seq_len(num_k)]),
      dmatrix1, dmatrix2
    )
    arreglar <- singletons.delete(table.groups.Pareto, population.pareto, num_k)
    population.pareto <- arreglar$poblacion
    population.R <- population.pareto

    if (isTRUE(local_search)) {
      population.R <- generate.path.relinking(population.pareto, num_k,
                                              dmatrix1, dmatrix2,
                                              number.objectives, cores)
      population.R <- generate.pareto.local.search(population.R, neighborhood,
                                                   num_k, num_objects, pop_size,
                                                   dmatrix1, dmatrix2,
                                                   number.objectives)
    }

    if (nrow(population.R) < pop_size) {
      population.P <- as.matrix(population.R[, seq_len(num_k)])
    } else {
      population.P <- as.matrix(population.R[seq_len(pop_size), seq_len(num_k)])
    }

    if (nrow(population.P) < pop_size) {
      population_random <- as.data.frame(
        t(sapply(seq_len(pop_size),
                 function(u) array(sample.int(num_objects, num_k, replace = FALSE))))
      )
      reparar <- singletons.repair(as.matrix(population_random),
                                   dmatrix1, dmatrix2, num_objects)
      population_random <- reparar$population
      if (nrow(population.P) > 1) {
        population.P <- population.P[!duplicated(population.P[, seq_len(num_k)]), ]
        population.P <- as.data.frame(population.P)
      }
      population.P <- rbind(population.P[, seq_len(num_k)],
                            population_random[(nrow(population.P) + 1):pop_size, ])
      rownames(population.P) <- seq_len(nrow(population.P))
      population.P <- as.matrix(population.P)
    }

    g <- g + 1
  }

  generate.results(num_k, dmatrix1, dmatrix2, pop_size, generation,
                   rat_cross, rat_muta, population.R)
}

#' @rdname moc.gapbk
#' @param ... Arguments passed to \code{\link{moc.gapbk}}.
#' @details
#' \code{moc.gabk} (note the single \code{p}) is a deprecated alias kept
#' for backward compatibility with versions 0.1.x. New code should call
#' \code{moc.gapbk} directly.
#' @export
moc.gabk <- function(...) {
  .Deprecated("moc.gapbk", package = "moc.gapbk",
              msg = paste0("'moc.gabk' is deprecated. ",
                           "Use 'moc.gapbk' instead. ",
                           "The alias will be removed in a future release."))
  moc.gapbk(...)
}
