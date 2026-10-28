#' @rdname efp2-methods
#' @method print efp2
#' @exportS3Method print efp2
print.efp2 <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
