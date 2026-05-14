# Internal helper: Path-Relinking (PR).
#
# Uses foreach::`%dopar%` registered through doParallel.
# The MPI backend (doMPI) has been removed for CRAN compliance because
# doMPI requires Rmpi, which depends on a system-level MPI installation
# and currently fails CRAN checks on macOS.

#' @importFrom foreach %dopar%
NULL

generate.path.relinking <- function(population_mejorar, num_k, dmatrix1,
                                    dmatrix2, number.objectives, cores) {

  i <- NULL  # silence R CMD check NOTE inside %dopar%

  lista_funciones <- c("generate.initial.population",
                       "generate.groups",
                       "verify.feasibility",
                       "singletons.delete",
                       "singletons.repair",
                       "generate.crossover.k.points",
                       "generate.mutation.random.controller",
                       "calculate.objective.functions",
                       "calculate.ranking.crowding",
                       "calculate.ranking",
                       "calculate.objectives.range",
                       "generate.results")

  lista_paquetes <- c("stats", "nsga2R",
                      "foreach", "parallel", "doParallel", "utils")

  if (nrow(population_mejorar) > 1) {

    if (nrow(population_mejorar) == 2) {
      pares_normal  <- t(c(1, 2))
      pares_inverso <- t(c(2, 1))
    } else {
      pares_normal  <- t(utils::combn(nrow(population_mejorar), 2))
      pares_inverso <- pares_normal[, c(2, 1)]
    }
    pares <- rbind(pares_normal, pares_inverso)

    soluciones_path <- data.frame()

    mopr <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(mopr)
    on.exit(parallel::stopCluster(mopr), add = TRUE)

    iteraciones_path <- foreach::foreach(
      i = seq_len(nrow(pares)),
      .export = lista_funciones,
      .packages = lista_paquetes,
      .combine = rbind
    ) %dopar% {

      invisible(rbind)
      soluciones_path <- data.frame()

      s_initial <- population_mejorar[pares[i, 1], seq_len(num_k)]
      s_guide   <- population_mejorar[pares[i, 2], seq_len(num_k)]

      # Cache the symmetric differences once per iteration of the while
      # (the original code recomputed setdiff three times per cycle).
      diff_si_sg <- setdiff(s_initial, s_guide)

      while (length(diff_si_sg) != 0) {
        elementos_eliminar <- colnames(diff_si_sg)
        elementos_agregar  <- colnames(setdiff(s_guide, s_initial))

        n_movs <- length(diff_si_sg)
        # Build n_movs intermediate solutions in one go by replicating
        # s_initial and overwriting one column per row. This avoids the
        # O(n_movs) chain of rbind() calls on a growing data frame.
        soluciones_intermedias <- s_initial[rep(1L, n_movs), , drop = FALSE]
        for (mov in seq_len(n_movs)) {
          soluciones_intermedias[mov, elementos_eliminar[mov]] <-
            s_guide[1L, elementos_agregar[mov]]
        }

        if (nrow(soluciones_intermedias) > 1) {
          tablagroupsIntermedias <- generate.groups(
            nrow(soluciones_intermedias),
            as.matrix(soluciones_intermedias[, seq_len(num_k)]),
            dmatrix1, dmatrix2
          )
          arreglar <- singletons.delete(tablagroupsIntermedias,
                                        soluciones_intermedias, num_k)
          soluciones_intermedias <- arreglar$poblacion
          tablagroupsIntermedias <- arreglar$groups

          if (length(arreglar$groups) > 0) {
            dominanciasIntermedias <- calculate.ranking.crowding(
              nrow(soluciones_intermedias),
              as.matrix(soluciones_intermedias[, seq_len(num_k)]),
              tablagroupsIntermedias, number.objectives, FALSE,
              dmatrix1, dmatrix2, num_k
            )
            dominanciasIntermedias <- as.data.frame(dominanciasIntermedias)

            soluciones_path <- rbind(soluciones_path,
                                     dominanciasIntermedias[1, seq_len(num_k)])
            s_initial <- dominanciasIntermedias[1, seq_len(num_k)]
          } else {
            s_initial <- s_guide
          }
        } else {
          tablagroupsIntermedias <- generate.groups(
            nrow(soluciones_intermedias),
            as.matrix(soluciones_intermedias[, seq_len(num_k)]),
            dmatrix1, dmatrix2
          )
          arreglar <- singletons.delete(tablagroupsIntermedias,
                                        soluciones_intermedias, num_k)
          soluciones_intermedias <- arreglar$poblacion
          tablagroupsIntermedias <- arreglar$groups

          if (length(arreglar$groups) > 0) {
            s_initial <- soluciones_intermedias[1, seq_len(num_k)]
          } else {
            s_initial <- s_guide
          }
        }
        # Refresh the cached symmetric difference for the loop condition.
        diff_si_sg <- setdiff(s_initial, s_guide)
      }

      soluciones_path
    }

    if (nrow(iteraciones_path) > 0) {
      soluciones_path <- as.matrix(iteraciones_path)
      soluciones_merge <- rbind(population_mejorar[, seq_len(num_k)],
                                soluciones_path[, seq_len(num_k)])
    } else {
      soluciones_merge <- population_mejorar[, seq_len(num_k)]
    }

    tablagroupsMerge <- generate.groups(
      nrow(soluciones_merge),
      as.matrix(soluciones_merge[, seq_len(num_k)]),
      dmatrix1, dmatrix2
    )
    arreglar <- singletons.delete(tablagroupsMerge,
                                  as.matrix(soluciones_merge[, seq_len(num_k)]),
                                  num_k)
    soluciones_merge <- arreglar$poblacion
    tablagroupsMerge <- arreglar$groups

    if (length(arreglar$groups) > 0) {
      dominanciasMerge <- calculate.ranking.crowding(
        nrow(soluciones_merge),
        as.matrix(soluciones_merge[, seq_len(num_k)]),
        tablagroupsMerge, number.objectives, FALSE,
        dmatrix1, dmatrix2, num_k
      )
      dominanciasMerge <- as.data.frame(dominanciasMerge)
      soluciones_no_dominadas <- subset(dominanciasMerge,
                                        dominanciasMerge$paretoranking == "1")
      soluciones_no_dominadas <- soluciones_no_dominadas[
        !duplicated(soluciones_no_dominadas[, seq_len(num_k)]),
      ]
      rownames(soluciones_no_dominadas) <- seq_len(nrow(soluciones_no_dominadas))
    } else {
      soluciones_no_dominadas <- population_mejorar
    }

    return(soluciones_no_dominadas)

  } else {
    warning("It is not possible to apply Path-Relinking (only one solution available).")
    rownames(population_mejorar) <- seq_len(nrow(population_mejorar))
    return(population_mejorar)
  }
}
