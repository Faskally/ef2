#' Conversion of an efp2 fit to a gam fit for plotting
#'
#' Converts an efp2 fit to a gam object, so that the plotting and prediction
#' tools in the \code{mgcv} package can be used. The coefficients and their
#' variances are those estimated by \code{efp2}, but the smoothing parameters
#' and effective degrees of freedom in the gam object are not.
#'
#' The conversion fits a model using \code{mgcv}, which can take a very
#' long time if the model has many random effects. Predictions do not need
#' the conversion, see \code{\link{predict.efp2}}.
#'
#'
#' @param object an efp2 fitted model
#' @return an object of class \code{gam}, see \code{\link[mgcv]{gamObject}}
#'
#' @examples
#' data(ef2_data)
#' ef2_data$sampleID <- factor(paste(ef2_data$date, ef2_data$siteID,
#'                                  ef2_data$lifestage, ef2_data$species))
#' m1 <- efp2(count ~ lifestage + s(doy, k = 3),
#'           data = ef2_data, pass = pass, id = sampleID)
#'
#' m1_gam <- as.gam2(m1)
#' plot(m1_gam)
#'
#' @importFrom mgcv gam
#' @importFrom stats binomial
#' @importFrom methods is
#'
#' @export
as.gam2 <- function(object) {
  stopifnot(is(object)[1] == "efp2")
  g1 <- efp_gam(object)

  # columns of the design matrix that were estimated
  whichkeep <- object$gam_cols
  redundant <- names(g1$coefficients[-1 * whichkeep])
  if (length(redundant)) {
    warning(
      "there is some redundancy in the model specification.",
      "This may be fine, but best to check."
    )
  }

  message("Note that, although predictions and plots will be correct\n the smoothing parameters and degrees of freedom in the gam object are not those estimated by efp2!")
  g1
}

# gam container holding the efp2 estimates, used for prediction and plotting
efp_gam <- function(object) {
  # make a gam container from the original setup. The smoothing parameters
  # are fixed so that this is quick: the gam fit itself is not used.
  G <- object $ Gsetup
  G$sp[] <- 1
  g1 <- gam(G = G)

  # columns of the design matrix that were estimated
  whichkeep <- object$gam_cols

  # assign coefficients to gam container
  g1$coefficients[] <- object$gam_coefficients

  # assign variance matrix to gam container
  g1$Vp[] <- 0
  diag(g1$Vp[]) <- 1e-5
  if (is.matrix(object$Vp)) {
    g1$Vp[whichkeep, whichkeep] <- object$Vp
  }

  # assign binomial family to gam container
  g1$family <- binomial()

  g1
}
