#' @importFrom stats pnorm qnorm
NULL

# Suppress R CMD CHECK notes for ggplot2 column names used as bare symbols
# inside aes() in manhattan_plot().
utils::globalVariables(c("gpos", "log10p", "colour_group"))
