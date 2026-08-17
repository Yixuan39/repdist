# G3 gate (method_notes.md §0, §3): weighted MMD must not depend on n_eff at fixed
# composition. The V-statistic does; the diagonal correction largely removes it;
# evenness matching removes it properly.

rdirichlet <- function(n, conc) {
  x <- matrix(rgamma(n * length(conc), conc), n, length(conc), byrow = TRUE)
  x / rowSums(x)
}

n_eff <- function(P) 1 / rowSums(P^2)  # effective support size; no longer exported

mmd2_v <- function(p, q, K) drop(p %*% K %*% p + q %*% K %*% q - 2 * p %*% K %*% q)

mmd2_u <- function(p, q, K) {
  dk <- diag(K)
  drop((p %*% K %*% p - sum(p^2 * dk)) / (1 - sum(p^2)) +
       (q %*% K %*% q - sum(q^2 * dk)) / (1 - sum(q^2)) - 2 * p %*% K %*% q)
}

alpha_transform <- function(p, alpha) { w <- p^alpha; w / sum(w) }

# n_eff is monotone DECREASING in alpha, so matching every sample DOWN to the
# minimum needs alpha >= 1 (sharpening). Sharpening only de-weights the
# low-abundance tail; matching upward toward uniform would inflate the influence
# of the least reliable measurements, so downward is the conservative direction.
match_n_eff <- function(p, target, hi = 64, iters = 80) {
  lo <- 1
  for (i in seq_len(iters)) {
    mid <- 0.5 * (lo + hi)
    if (n_eff(rbind(alpha_transform(p, mid))) > target) lo <- mid else hi <- mid
  }
  0.5 * (lo + hi)
}

set.seed(0)
U <- 400L; d <- 16L
Zc <- matrix(rnorm(U * d), U, d)
Zc <- Zc / sqrt(rowSums(Zc^2))
Cc <- as.matrix(dist(Zc))
Kc <- rbf_kernel(Cc)
k_bar <- (sum(Kc) - sum(diag(Kc))) / (U * (U - 1))

test_that("n_eff is monotone in the sharpening exponent", {
  p <- drop(rdirichlet(1, rep(0.3, U)))
  ne <- vapply(c(0.1, 0.25, 0.5, 1, 2, 4),
               function(a) n_eff(rbind(alpha_transform(p, a))), numeric(1))
  expect_true(all(diff(ne) < 0))
})

test_that("the V-statistic tracks the §3 excess-concentration prediction", {
  # Fixed shared universe (the metaproteomic case: the metagenome fixes the support),
  # so the driver is excess concentration above uniform, not 1/n_eff itself.
  kap <- function(p) sum(p^2) - 1 / U

  res <- t(vapply(c(0.05, 0.2, 1, 5), function(conc) {
    v <- u <- pred <- numeric(60)
    for (i in 1:60) {
      a <- drop(rdirichlet(1, rep(conc, U)))
      b <- drop(rdirichlet(1, rep(conc, U)))
      v[i] <- mmd2_v(a, b, Kc)
      u[i] <- mmd2_u(a, b, Kc)
      pred[i] <- (1 - k_bar) * (kap(a) + kap(b))
    }
    c(v = mean(v), u = mean(u), pred = mean(pred))
  }, numeric(3)))

  expect_lt(max(abs(res[, "v"] - res[, "pred"]) / res[, "pred"]), 0.15)
  # Contamination is real: concentrated samples sit further apart than even ones.
  expect_gt(res[1, "v"], 5 * res[nrow(res), "v"])

  # The U-statistic reduces but does not eliminate the dependence. The residual
  # offset -2(1-k_bar)/U is a fixed-universe artifact: constant across pairs, so
  # benign for PCoA and invariant under permutation.
  scale <- max(res[, "v"])
  v_spread <- diff(range(res[, "v"])) / scale
  u_spread <- diff(range(res[, "u"])) / scale
  expect_gt(v_spread / u_spread, 10)
  offset <- -2 * (1 - k_bar) / U
  expect_lt(abs(median(res[, "u"]) - offset), 0.2 * abs(offset))
})

test_that("evenness matching equalizes n_eff", {
  ps <- lapply(c(0.05, 0.2, 1, 5), function(c) drop(rdirichlet(1, rep(c, U))))
  target <- min(vapply(ps, function(x) n_eff(rbind(x)), numeric(1)))
  matched <- lapply(ps, function(x) alpha_transform(x, match_n_eff(x, target)))
  got <- vapply(matched, function(x) n_eff(rbind(x)), numeric(1))
  expect_lt(diff(range(got)), 0.01 * target)
})

test_that("the weighted MMD satisfies the §4 metric claims", {
  ps <- lapply(c(0.05, 0.2), function(c) drop(rdirichlet(1, rep(c, U))))
  expect_lt(abs(mmd2_v(ps[[1]], ps[[1]], Kc)), 1e-12)
  expect_equal(mmd2_v(ps[[1]], ps[[2]], Kc), mmd2_v(ps[[2]], ps[[1]], Kc))
  P <- do.call(rbind, ps)
  expect_gt(min(eigen(P %*% Kc %*% t(P), symmetric = TRUE, only.values = TRUE)$values),
            -1e-10)
})
