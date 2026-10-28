# Set up the design matrices of a model formula, using mgcv to construct
# the terms. Penalised terms (random effects and penalised smooths) are split
# into fixed effects, for the part of the term that is not penalised, and
# random effects, see smooth_setup.
#   formula: a gam formula
#   data: a data.frame
#   sparse: passed to smooth_setup
# Returns a list with the mgcv setup Gsetup, its design matrix G, the fixed
# and random effects design matrices X and Z, the penalised blocks, the
# labels of their standard deviations, and the matrix to_gam which converts
# the parameters (fixed then random effects) to the coefficients of G.
model_setup <- function(formula, data, sparse = TRUE) {
  Gsetup <- gam(formula, data = data, fit = FALSE, drop.unused.levels = FALSE)
  G <- Gsetup$X
  colnames(G) <- Gsetup$term.names
  if (nrow(G) != nrow(data)) stop("I think there are NAs in your data, please check and remove them before rerunning.")

  # identify random effects: these are the s(..., bs = "re") terms
  is_re <- vapply(Gsetup$smooth, inherits, TRUE, what = "random.effect")
  re_cols <- lapply(Gsetup$smooth[is_re], function(x) x$first.para:x$last.para)
  re_cols <- as.integer(unlist(re_cols))

  # set up the penalised terms: random effects and penalised smooths. Each
  # is split into fixed effects, for the part of the term that is not
  # penalised, and random effects
  blocks <- smooth_setup(Gsetup, formula, sparse = sparse)
  pen_cols <- as.integer(unlist(lapply(blocks, function(x) x$cols)))
  fx_cols <- setdiff(seq_len(ncol(G)), pen_cols)
  nsigma <- vapply(blocks, function(x) x$npen, 0)
  re_labels <-
    unlist(lapply(blocks, function(x) {
      if (x$npen == 1) x$label else paste0(x$label, seq_len(x$npen))
    }))
  if (is.null(re_labels)) re_labels <- character(0)

  # the fixed effects design matrix
  X <- do.call(cbind, c(list(G[, fx_cols, drop = FALSE]), lapply(blocks, function(x) x$Xf)))
  Z <- do.call(cbind, c(list(Matrix(0, nrow(G), 0)), lapply(blocks, function(x) x$Zu)))

  # remove redundant / collinear parameters (fixed effects only)
  qr.X <- qr(X)
  rank.deficient <- qr.X $ pivot[abs(diag(qr.X $ qr)) < 1e-7]
  if (length(rank.deficient)) {
    droppar <- paste(colnames(X)[rank.deficient], collapse = "\n\t")
    warning("*** Model has ", length(rank.deficient)," too many parameter(s)!!\n    i will remove the redundant ones:\n\t", droppar, call. = FALSE)
  }
  fx_keep <- setdiff(seq_len(ncol(X)), rank.deficient)
  fx_cols <- setdiff(fx_cols, fx_cols[rank.deficient[rank.deficient <= length(fx_cols)]])

  # the matrix which converts the model parameters (fixed then random
  # effects) to the coefficients of the columns of G
  nfx <- ncol(G) - length(pen_cols)
  nX <- ncol(X)
  triplets <- function(M, rows, cols) {
    M <- as.matrix(M)
    nz <- which(M != 0, arr.ind = TRUE)
    cbind(rows[nz[, 1]], cols[nz[, 2]], M[nz])
  }
  to_gam <- list(cbind(setdiff(seq_len(ncol(G)), pen_cols), seq_len(nfx), rep(1, nfx)))
  at_g <- nfx
  at_u <- nX
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    ipar <- c(at_g + seq_len(ncol(block$Xf)), at_u + seq_len(ncol(block$Zu)))
    to_gam <- c(to_gam, list(triplets(block$TG, block$cols, ipar)))
    if (!is.null(block$Tint)) {
      to_gam <- c(to_gam, list(triplets(t(block$Tint), which(colnames(G) == "(Intercept)"), ipar)))
    }
    # where the random effects and standard deviations are
    blocks[[i]]$idx <- at_u - nX + seq_len(ncol(block$Zu))
    blocks[[i]]$sig <- sum(nsigma[seq_len(i - 1)]) + seq_len(block$npen)
    at_g <- at_g + ncol(block$Xf)
    at_u <- at_u + ncol(block$Zu)
  }
  to_gam <- do.call(rbind, to_gam)
  to_gam <-
    sparseMatrix(
      i = to_gam[, 1], j = to_gam[, 2], x = to_gam[, 3],
      dims = c(ncol(G), nX + ncol(Z))
    )
  to_gam <- to_gam[, c(fx_keep, nX + seq_len(ncol(Z))), drop = FALSE]
  gam_cols <- which(rowSums(abs(to_gam)) > 0)
  Gfit <- G[, gam_cols, drop = FALSE]

  list(
    Gsetup = Gsetup, G = G, X = X, Z = Z, blocks = blocks,
    re_labels = re_labels, re_cols = re_cols, fx_cols = fx_cols,
    fx_keep = fx_keep, to_gam = to_gam, gam_cols = gam_cols, Gfit = Gfit
  )
}

# the parts of a model set up by model_setup that are needed by the
# likelihood
model_data <- function(ms, offset) {
  list(
    offset = offset,
    X = as(as(ms$X[, ms$fx_keep, drop = FALSE], "generalMatrix"), "TsparseMatrix"),
    Z = as(as(ms$Z, "generalMatrix"), "TsparseMatrix"),
    blocks = lapply(ms$blocks, function(x) x[c("idx", "sig", "S")])
  )
}
