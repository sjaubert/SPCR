## sim_abc.R : simulation de la section 6.3 du document ANOVA à plusieurs facteurs (environ 1 min)
## Simulation : négliger ABC (modèle fixe 3x3x3, r = 1) quand ABC existe vraiment
set.seed(2026)
lev <- expand.grid(A=factor(1:3), B=factor(1:3), C=factor(1:3))
## motif ABC fixe : projection d'un tirage sur l'espace de l'interaction triple
z <- rnorm(27)
abc0 <- residuals(lm(z ~ (A+B+C)^2, data = cbind(lev, z = z)))   # sommes nulles sur chaque indice
abc0 <- abc0 / sqrt(sum(abc0^2)/8)     # sum(abc^2)/ddl_ABC = 1
## vérification des contraintes
stopifnot(max(abs(tapply(abc0, list(lev$A, lev$B), sum))) < 1e-12)
sim <- function(k, effA = 0, nsim = 4000) {
  rej <- 0
  for (s in 1:nsim) {
    y <- effA * (as.integer(lev$A) - 2) + k * abc0 + rnorm(27)
    p <- summary(aov(y ~ (A+B+C)^2, data = cbind(lev, y = y)))[[1]][1, "Pr(>F)"]
    rej <- rej + (p < 0.05)
  }
  rej / nsim
}
for (k in c(0, 1, 2)) cat(sprintf("k=%d  E(CM_ABC)=%.0f  rejet A nul %.3f  puissance effA=1 %.3f\n",
                                  k, 1 + k^2, sim(k), sim(k, 1)))
