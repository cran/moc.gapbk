## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----basic--------------------------------------------------------------------
library(moc.gapbk)

set.seed(2025)

# Toy data: 50 objects (e.g. genes) described by 20 features (e.g. samples).
x <- matrix(stats::runif(50 * 20, min = -5, max = 10),
            nrow = 50, ncol = 20)

# Two distance matrices over the same set of objects.
# Here we use amap if available (correlation distance is biologically
# common), and fall back to base R otherwise so the vignette knits
# under any configuration.
if (requireNamespace("amap", quietly = TRUE)) {
  d1 <- as.matrix(amap::Dist(x, method = "euclidean"))
  d2 <- as.matrix(amap::Dist(x, method = "correlation"))
} else {
  d1 <- as.matrix(stats::dist(x, method = "euclidean"))
  d2 <- as.matrix(stats::dist(x, method = "manhattan"))
}

res <- moc.gapbk(dmatrix1 = d1,
                 dmatrix2 = d2,
                 num_k = 3,
                 generation = 5,
                 pop_size = 6)

## ----population---------------------------------------------------------------
head(res$population)

## ----matrix-solutions---------------------------------------------------------
head(res$matrix.solutions)

## ----clustering-vec-----------------------------------------------------------
str(res$clustering[[1]])
table(res$clustering[[1]])

## ----local-search, eval = FALSE-----------------------------------------------
# res_full <- moc.gapbk(d1, d2,
#                       num_k = 3,
#                       generation = 10,
#                       pop_size = 10,
#                       local_search = TRUE,
#                       cores = 2)

