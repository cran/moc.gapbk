# Internal helper: assemble final results returned by moc.gapbk().

generate.results <- function(num_k, dmatrix1, dmatrix2, pop_size,
                             generation, rat_cross, rat_muta, population) {

  d <- as.data.frame(population[seq_len(pop_size), ])
  D <- d[order(d$paretoranking, -d$crowding), ]
  paretos <- subset(D, D$paretoranking == "1")
  nondom <- paretos

  population.PARETO <- as.matrix(nondom[, seq_len(num_k)])
  table.groups.PARETO <- generate.groups(nrow(population.PARETO),
                                         population.PARETO, dmatrix1, dmatrix2)

  pareto_solutions <- as.data.frame(table.groups.PARETO)
  colnames(pareto_solutions) <- seq_len(length(table.groups.PARETO))

  vectors.partition.list <- list()
  for (i in seq_len(ncol(pareto_solutions))) {
    vector.partition <- unlist(pareto_solutions[, i])
    names(vector.partition) <- rownames(pareto_solutions)
    vectors.partition.list[[i]] <- vector.partition
  }

  list(population = nondom,
       matrix.solutions = pareto_solutions,
       clustering = vectors.partition.list)
}
