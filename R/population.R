# Internal helpers: population construction and feasibility repair.
# Not exported.

generate.initial.population <- function(num_objects, num_k, pop_size) {
  population.P <- t(sapply(seq_len(pop_size),
                           function(x) sample.int(num_objects, num_k, replace = FALSE)))
  as.matrix(population.P)
}

verify.feasibility <- function(population, par_destino, num_objects) {
  # For each individual, replace duplicated medoids at random
  for (p in seq_len(nrow(population))) {
    while (length(which(duplicated(population[p, ]))) > 0) {
      id_rep <- which(duplicated(population[p, ]))
      for (i in seq_along(id_rep)) {
        population[p, id_rep[i]] <- sample.int(num_objects, 1, replace = FALSE)
      }
    }
  }
  population
}
