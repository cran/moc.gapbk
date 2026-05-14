# Internal helpers: cluster formation and singleton handling.
#
# Vectorised since 0.3.0: the original nested for-loop in
# generate.groups() has been replaced by an O(n*k) matrix lookup, and
# singletons.delete()/singletons.repair() use vectorised counts.

# Assign each object to its nearest medoid for every individual in the
# population.
#
# medoids: matrix (pop_size x k) of medoid indices.
# Returns a list of integer vectors of length n (one per individual).
generate.groups <- function(pop_size, medoids, dmatrix1, dmatrix2) {
  if (is.null(dim(medoids))) {
    medoids <- matrix(medoids, nrow = 1)
  }
  n <- nrow(dmatrix1)
  groups <- vector("list", pop_size)
  rn <- rownames(dmatrix1)

  for (p in seq_len(pop_size)) {
    # dmatrix1[, medoids[p, ]] is an (n x k) submatrix; the cluster of
    # object i is which.min over its row. max.col(-X) is the vectorised
    # rowwise which.min and runs in C (orders of magnitude faster than
    # the original double for loop).
    sub <- dmatrix1[, medoids[p, ], drop = FALSE]
    formacion <- max.col(-sub, ties.method = "first")
    names(formacion) <- rn
    groups[[p]] <- formacion
  }
  groups
}

# Drop individuals that contain at least one singleton cluster.
# Vectorised count via tabulate() rather than which()+length().
singletons.delete <- function(groups_formados, poblacion, num_k) {
  if (length(groups_formados) == 0L) {
    return(list(groups = groups_formados, poblacion = poblacion))
  }

  has_singleton <- vapply(groups_formados, function(g) {
    tab <- tabulate(g, nbins = num_k)
    any(tab < 2L)
  }, logical(1))

  if (any(has_singleton)) {
    keep <- !has_singleton
    groups_formados <- groups_formados[keep]
    if (is.matrix(poblacion) || is.data.frame(poblacion)) {
      poblacion <- poblacion[keep, , drop = FALSE]
    } else {
      poblacion <- poblacion[keep]
    }
  }

  list(groups = groups_formados, poblacion = poblacion)
}

# Repair singletons by sampling a fresh medoid from the unused pool.
# The original implementation reset the column index to 1 every time a
# replacement happened, leading to O(n^2) worst case retries. The new
# version retries only at the same column and bounds the number of
# attempts.
singletons.repair <- function(population.repair, dmatrix1, dmatrix2, num_objects) {
  if (is.null(dim(population.repair))) {
    population.repair <- matrix(population.repair, nrow = 1)
  }

  table.repair.groups <- generate.groups(nrow(population.repair),
                                         population.repair, dmatrix1, dmatrix2)

  max_attempts <- num_objects  # safety bound
  num_k <- ncol(population.repair)

  for (individuo in seq_len(nrow(population.repair))) {
    for (columna in seq_len(num_k)) {
      attempts <- 0L
      # current cluster sizes for this individual
      sizes <- tabulate(table.repair.groups[[individuo]], nbins = num_k)
      while (sizes[columna] < 3L && attempts < max_attempts) {
        genes_cromosoma <- population.repair[individuo, ]
        posibles <- setdiff(seq_len(num_objects), genes_cromosoma)
        if (length(posibles) == 0L) break
        reemplazo <- if (length(posibles) == 1L) posibles
                     else sample(posibles, 1L)
        population.repair[individuo, columna] <- reemplazo
        new_groups <- generate.groups(1L,
                                      population.repair[individuo, , drop = FALSE],
                                      dmatrix1, dmatrix2)
        table.repair.groups[[individuo]] <- new_groups[[1L]]
        sizes <- tabulate(table.repair.groups[[individuo]], nbins = num_k)
        attempts <- attempts + 1L
      }
    }
  }

  list(groups = table.repair.groups, population = population.repair)
}
