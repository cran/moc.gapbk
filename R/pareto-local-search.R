# Internal helper: Pareto Local Search (PLS).
#
# Optimised since 0.3.0: row-matching that used to be
# `apply(A, 1, function(x) all(x == B))` is now a single
# `rowSums(A != B) == 0L`, which runs in vectorised C and avoids
# allocating a closure per row. Combined with the vectorised
# generate.groups() and calculate.objective.functions() this changes
# the inner-loop cost from O(m*k) per probe (R-level) to O(m*k) (C-level).

row_equals <- function(A, B) {
  # Returns logical indicating which rows of A equal the single row B.
  if (is.vector(A)) A <- matrix(A, nrow = 1L)
  if (is.vector(B)) B <- matrix(B, nrow = 1L)
  B_vec <- as.numeric(B[1L, ])
  # rowSums in C is much faster than apply(_, 1, all).
  rowSums(A != matrix(B_vec, nrow = nrow(A), ncol = length(B_vec),
                      byrow = TRUE)) == 0L
}

generate.pareto.local.search <- function(population_pareto, neighborhood,
                                         num_k, num_objects, pop_size,
                                         dmatrix1, dmatrix2, number.objectives) {

  if (nrow(population_pareto) > pop_size) {
    poblacion_A0 <- population_pareto[seq_len(pop_size), seq_len(num_k)]
  } else {
    poblacion_A0 <- population_pareto[, seq_len(num_k)]
  }

  poblacion_A0 <- poblacion_A0[!duplicated(poblacion_A0[, seq_len(num_k)]), ]
  poblacion_A0 <- cbind(explored = rep(0, nrow(poblacion_A0)), poblacion_A0)
  names(poblacion_A0) <- paste("V", seq_len(num_k + 1), sep = "")
  poblacion_A <- poblacion_A0

  num_neighborhood <- ceiling(neighborhood * num_objects)

  if (nrow(poblacion_A0) > 1) {

    while (nrow(poblacion_A0) > 1) {
      posicion_alterar <- sample.int(num_k, 1)
      s <- sample.int(nrow(poblacion_A0), 1)
      solucion_S <- as.data.frame(poblacion_A0[s, ])

      lista_medoides <- as.vector(as.matrix(poblacion_A0[, 2:(num_k + 1)]))
      lista_medoides_unicos <- unique(lista_medoides)
      lista_medoides_total <- seq_len(num_objects)
      lista_medoides_disponibles <- setdiff(lista_medoides_total,
                                            lista_medoides_unicos)

      for (i in seq_len(num_neighborhood)) {
        if (length(lista_medoides_disponibles) < 1) break
        nuevo_medoide <- lista_medoides_disponibles[
          sample.int(length(lista_medoides_disponibles), 1, replace = FALSE)
        ]
        lista_medoides_disponibles <- setdiff(lista_medoides_disponibles,
                                              nuevo_medoide)

        if (length(lista_medoides_disponibles) > 0) {
          solucion_S_prima <- solucion_S
          solucion_S_prima[1, (posicion_alterar + 1)] <- nuevo_medoide

          tablagroupsS_prima <- generate.groups(
            nrow(solucion_S_prima),
            as.matrix(solucion_S_prima[, 2:(num_k + 1)]),
            dmatrix1, dmatrix2
          )
          arreglar <- singletons.delete(tablagroupsS_prima, solucion_S_prima, num_k)

          if (length(arreglar$groups) > 0) {
            poblacion_unida <- rbind(solucion_S_prima, poblacion_A)
            poblacion_unida <- poblacion_unida[
              !duplicated(poblacion_unida[, 2:(num_k + 1)]),
            ]

            table.groups.PU <- generate.groups(
              nrow(poblacion_unida),
              as.matrix(poblacion_unida[, 2:(num_k + 1)]),
              dmatrix1, dmatrix2
            )
            arreglar <- singletons.delete(table.groups.PU, poblacion_unida, num_k)
            poblacion_unida <- arreglar$poblacion
            table.groups.PU <- arreglar$groups

            dominanciasPU <- calculate.ranking.crowding(
              nrow(poblacion_unida),
              as.matrix(poblacion_unida[, seq_len(num_k + 1)]),
              table.groups.PU, number.objectives, TRUE,
              dmatrix1, dmatrix2, num_k
            )
            dominanciasPU <- as.data.frame(dominanciasPU)

            calidad_solucion_S_prima <- merge(dominanciasPU, solucion_S_prima)

            if (is.na(calidad_solucion_S_prima$crowding)) {
              calidad_solucion_S_prima$crowding <- 0
            }
            if (is.na(dominanciasPU$crowding[1])) dominanciasPU$crowding[1] <- 0

            if (calidad_solucion_S_prima$paretoranking == "1" &&
                calidad_solucion_S_prima$crowding >= dominanciasPU$crowding[1]) {

              A <- dominanciasPU[, 2:(num_k + 1)]
              B <- solucion_S_prima[, 2:(num_k + 1)]

              # Vectorised replacement of apply(A, 1, function(x) all(x == B)).
              fila_prima <- which(row_equals(as.matrix(A), as.matrix(B)))
              if (length(fila_prima) > 0L) {
                dominanciasPU[fila_prima, 1L] <- 0
              }

              solo_las_no_dominadas <- as.data.frame(
                subset(dominanciasPU, dominanciasPU$paretoranking == "1")
              )
              poblacion_A <- solo_las_no_dominadas[, seq_len(num_k + 1)]
            }
          }
        }
      }

      solucion_S <- as.data.frame(solucion_S)
      J <- poblacion_A[, 2:(num_k + 1)]
      P <- solucion_S[, 2:(num_k + 1)]

      fila_solucion_s <- which(row_equals(as.matrix(J), as.matrix(P)))
      if (length(fila_solucion_s) > 0) {
        poblacion_A[fila_solucion_s, 1L] <- 1
      }

      poblacion_A <- as.data.frame(poblacion_A)
      poblacion_A0 <- as.data.frame(subset(poblacion_A, poblacion_A[, 1] == "0"))
    }

    poblacion_A <- poblacion_A[, 2:(num_k + 1)]
    colnames(poblacion_A) <- paste("V", seq_len(num_k), sep = "")
    poblacion_A <- poblacion_A[!duplicated(poblacion_A[, seq_len(num_k)]), ]
    poblacion_A <- as.data.frame(poblacion_A)

    tablagroupsA <- generate.groups(
      nrow(poblacion_A[, seq_len(num_k)]),
      as.matrix(poblacion_A[, seq_len(num_k)]),
      dmatrix1, dmatrix2
    )
    poblacion_A <- calculate.ranking.crowding(
      nrow(poblacion_A[, seq_len(num_k)]),
      as.matrix(poblacion_A[, seq_len(num_k)]),
      tablagroupsA, number.objectives, FALSE,
      dmatrix1, dmatrix2, num_k
    )
    poblacion_A <- as.data.frame(poblacion_A)

    return(as.data.frame(poblacion_A[, ]))

  } else {
    rownames(population_pareto) <- seq_len(nrow(population_pareto))
    return(population_pareto[1, ])
  }
}
