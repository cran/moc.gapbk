#' moc.gapbk: Multi-Objective Clustering Guided by a-Priori Biological Knowledge
#'
#' The \pkg{moc.gapbk} package implements the MOC-GaPBK algorithm proposed
#' by Parraga-Alava and others (2018). It combines NSGA-II with
#' Path-Relinking and Pareto Local Search to discover clustering
#' solutions that are good with respect to two objective functions
#' simultaneously, typically defined from two distance matrices: one
#' over the data itself and one encoding a-priori biological knowledge.
#'
#' The main user-facing function is \code{\link{moc.gapbk}}. The legacy
#' name \code{\link{moc.gabk}} is preserved as a deprecated alias for
#' backward compatibility.
#'
#' @references
#' J. Parraga-Alava, M. Dorn, M. Inostroza-Ponta (2018).
#' \emph{A multi-objective gene clustering algorithm guided by apriori
#' biological knowledge with intensification and diversification
#' strategies}. BioData Mining. 11(1) 1-16.
#' \doi{10.1186/s13040-018-0178-4}.
#'
#' @keywords internal
"_PACKAGE"

# Avoid R CMD check NOTE "no visible binding for global variable 'i'"
# triggered by the foreach %dopar% loop in generate.path.relinking().
utils::globalVariables(c("i"))
