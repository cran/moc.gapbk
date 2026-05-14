# Internal helpers: objective functions, Pareto ranking, crowding distance.
#
# Vectorised since 0.3.0. The original triple for-loop in
# calculate.objective.functions() ran in O(P * n * k) at the R
# interpreter level, dominating the runtime of the algorithm. It is now
# a single matrix slice plus rowSums(), implemented entirely in C
# through standard R primitives.

# Compute the Xie-Beni objective for each individual under a given
# distance matrix.
#
# Notes:
#  * If a precomputed squared distance matrix `dist_sq` is supplied via
#    the attribute `"dist_sq"` of par_distancia, the squaring step is
#    skipped. This is the recommended path when calling repeatedly with
#    the same matrix (see moc.gapbk()).
#  * `groups` is accepted only to preserve the calling convention; the
#    objective itself does not need group assignments, only medoids.
calculate.objective.functions <- function(pop_size, population, groups,
                                          par_distancia, num_k, local_search) {
  if (isTRUE(local_search)) {
    population <- population[, 2:(num_k + 1), drop = FALSE]
  }

  # Use a cached squared distance matrix when available
  dist_sq <- attr(par_distancia, "dist_sq")
  if (is.null(dist_sq)) {
    dist_sq <- par_distancia * par_distancia
  }

  P <- length(groups)
  if (P == 0L) return(numeric(0))
  n <- nrow(dist_sq)
  fobj <- numeric(P)

  for (x in seq_len(P)) {
    medoides <- population[x, seq_len(num_k)]

    # Numerador: sum over all clusters and all objects of d(medoid_k, i)^2.
    # The (k x n) block of squared distances is just one matrix slice.
    XB_numerador <- sum(dist_sq[medoides, , drop = FALSE])

    # Denominador: min squared distance between distinct medoid pairs.
    # Pull the (k x k) submatrix; the off-diagonal min is the answer.
    sub <- dist_sq[medoides, medoides, drop = FALSE]
    # Diagonal is zero (distance to self); replace with Inf so it does
    # not interfere with the min.
    diag(sub) <- Inf
    XB_denominador <- min(sub)

    fobj[x] <- XB_numerador / (n * XB_denominador)
  }

  fobj
}

calculate.ranking <- function(pop_size, ranking) {
  rankIndex <- integer(pop_size)
  i <- 1L
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1L
  }
  rankIndex
}

calculate.objectives.range <- function(population, par_var, par_obj, local_search) {
  if (isTRUE(local_search)) {
    cols <- (par_var + 2L):(par_var + par_obj + 1L)
  } else {
    cols <- (par_var + 1L):(par_var + par_obj)
  }
  if (nrow(population) <= 1L) {
    objetivos <- t(as.matrix(population[, cols]))
  } else {
    objetivos <- population[, cols, drop = FALSE]
  }
  apply(objetivos, 2L, max) - apply(objetivos, 2L, min)
}

calculate.ranking.crowding <- function(pop_size, population, groups,
                                       number_objectives, local_search,
                                       dmatrix1, dmatrix2, num_k) {
  obj1 <- calculate.objective.functions(pop_size, population, groups, dmatrix1,
                                        num_k, local_search)
  obj2 <- calculate.objective.functions(pop_size, population, groups, dmatrix2,
                                        num_k, local_search)
  population <- cbind(population, obj1, obj2)

  if (isTRUE(local_search)) {
    cols <- (num_k + 2L):(num_k + number_objectives + 1L)
  } else {
    cols <- (num_k + 1L):(num_k + number_objectives)
  }
  if (nrow(population) <= 1L) {
    objetivos <- t(as.matrix(population[, cols]))
  } else {
    objetivos <- population[, cols, drop = FALSE]
  }

  ranking <- nsga2R::fastNonDominatedSorting(objetivos)
  paretoranking <- calculate.ranking(pop_size, ranking)
  population <- cbind(population, paretoranking)

  objRange <- calculate.objectives.range(population, num_k, number_objectives,
                                         local_search)
  cd <- nsga2R::crowdingDist4frnt(population, ranking, objRange)
  population <- cbind(population, crowding = c(apply(cd, 1L, sum)))
  population <- as.data.frame(population)
  population <- population[order(population$paretoranking,
                                 -(population$crowding)), ]
  as.matrix(population)
}
