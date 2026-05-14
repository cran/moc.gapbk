# Internal helpers: genetic operators (crossover and mutation).

generate.crossover.k.points <- function(pop_size, num_k, population.from,
                                        population.to, rat_cross) {
  p <- 0.50

  for (i in seq(1, pop_size, by = 2)) {
    parejas <- sample.int(pop_size, 2, replace = FALSE)
    cruce <- stats::runif(1, 0, 1)

    if (cruce < rat_cross) {
      for (j in seq_len(num_k)) {
        aleatorio <- stats::runif(1, 0, 1)
        if (aleatorio <= p) {
          population.to[i,       j] <- as.matrix(population.from[parejas[1], j])
          population.to[(i + 1), j] <- as.matrix(population.from[parejas[2], j])
        } else {
          population.to[i,       j] <- as.matrix(population.from[parejas[2], j])
          population.to[(i + 1), j] <- as.matrix(population.from[parejas[1], j])
        }
      }
    } else {
      population.to[i,       seq_len(num_k)] <- as.matrix(population.from[parejas[1], seq_len(num_k)])
      population.to[(i + 1), seq_len(num_k)] <- as.matrix(population.from[parejas[2], seq_len(num_k)])
    }
  }
  population.to
}

generate.mutation.random.controller <- function(pop_size, num_k, population,
                                                rat_muta, num_objects) {
  for (p in seq_len(pop_size)) {
    muta <- stats::runif(1, 0, 1)
    if (muta < rat_muta) {
      posicion_mutar <- sample.int(num_k, 1, replace = FALSE)
      medoide <- sample.int(num_objects, 1, replace = FALSE)
      population[p, posicion_mutar] <- medoide
    }
  }
  population
}
