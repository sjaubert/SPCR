## rr_etude.R : calculs du document « Gage R&R et ANOVA »
## S. Jaubert, novembre 2026. R >= 4.0, paquet lme4 (contrôle REML).
## À lancer depuis le dossier ANOVA_et_RR du dépôt SPCR.
source("rr_functions.R")
suppressMessages(library(lme4))
rr <- function(df, IT = NA, k = 6, alpha_int = 0.05, label = "") {
  df$P <- factor(df$P); df$O <- factor(df$O)
  p <- nlevels(df$P); o <- nlevels(df$O); r <- nrow(df) / (p * o)
  s <- summary(aov(y ~ P * O, data = df))[[1]]
  SP <- s[1, 3]; SO <- s[2, 3]; SPO <- s[3, 3]; SE <- s[4, 3]
  Fint <- SPO / SE; pint <- pf(Fint, (p-1)*(o-1), p*o*(r-1), lower.tail = FALSE)
  out <- list(label = label, p = p, o = o, r = r, tab = s, Fint = Fint, pint = pint)
  if (pint < alpha_int) {
    vE <- SE; vPO <- max(0, (SPO - SE) / r); vO <- max(0, (SO - SPO) / (p * r)); vP <- max(0, (SP - SPO) / (o * r))
    ci <- mls_interaction(SP, SO, SPO, SE, p, o, r); modele <- "avec interaction"
    FP <- SP / SPO; FO <- SO / SPO
  } else {
    SEp <- (s[3, 2] + s[4, 2]) / (s[3, 1] + s[4, 1]); vE <- SEp; vPO <- 0
    vO <- max(0, (SO - SEp) / (p * r)); vP <- max(0, (SP - SEp) / (o * r))
    ci <- mls_additive(SP, SO, SEp, p, o, r); modele <- "sans interaction (regroupée)"; SE <- SEp
    FP <- SP / SEp; FO <- SO / SEp
  }
  vR <- vO + vPO; vM <- vE + vR; vT <- vM + vP
  ndc <- sqrt(2) * sqrt(vP / vM)
  crit <- c(sE = sqrt(vE), sR = sqrt(vR), sM = sqrt(vM), sP = sqrt(vP), sT = sqrt(vT),
            contribM = 100 * vM / vT, studyM = 100 * sqrt(vM / vT),
            tolM = if (!is.na(IT)) 100 * k * sqrt(vM) / IT else NA, ndc = ndc, ndc_trunc = floor(ndc))
  ## REML (lme4) pour comparaison
  fm <- if (pint < alpha_int) lmer(y ~ 1 + (1|P) + (1|O) + (1|P:O), data = df) else lmer(y ~ 1 + (1|P) + (1|O), data = df)
  vc <- as.data.frame(VarCorr(fm))[, c("grp", "vcov")]
  c(out, list(modele = modele, vE = vE, vPO = vPO, vO = vO, vP = vP, crit = crit, ci = ci, reml = vc,
              FP = FP, FO = FO))
}
show <- function(x) { cat("=====", x$label, "p o r =", x$p, x$o, x$r, "\n"); print(x$tab); cat("F interaction", x$Fint, "p", x$pint, "->", x$modele, "\n")
  cat("composantes : E", x$vE, " PO", x$vPO, " O", x$vO, " P", x$vP, "\n"); print(round(x$crit, 5)); print(x$ci); print(x$reml); cat("FP", x$FP, "FO", x$FO, "\n") }

set.seed(20261104)
p <- 10; o <- 3; r <- 3
X <- 11.925 + rnorm(p, 0, 0.020)                 # vraies valeurs des butées
bO_mic <- rnorm(o, 0, 0.0012); bPO_mic <- matrix(rnorm(p*o, 0, 0.0006), p, o)
bO_pied <- rnorm(o, 0, 0.004);  bPO_pied <- matrix(rnorm(p*o, 0, 0.002), p, o)
grid <- expand.grid(k = 1:r, O = 1:o, P = 1:p)
mic <- with(grid, round(X[P] + bO_mic[O] + bPO_mic[cbind(P, O)] + rnorm(nrow(grid), 0, 0.0015), 3))
pied <- with(grid, round((X[P] + bO_pied[O] + bPO_pied[cbind(P, O)] + rnorm(nrow(grid), 0, 0.006)) / 0.02) * 0.02)
tp_mic <- data.frame(P = grid$P, O = grid$O, k = grid$k, y = mic)
tp_pied <- data.frame(P = grid$P, O = grid$O, k = grid$k, y = round(pied, 2))
write.csv2(tp_mic, "tp_micrometre.csv", row.names = FALSE); write.csv2(tp_pied, "tp_pied.csv", row.names = FALSE)
cat("vraies valeurs X :", round(X, 4), " sd(X)", sd(X), "\n")
A <- rr(tp_mic, IT = 0.15, label = "TP micrometre"); show(A)
B <- rr(tp_pied, IT = 0.15, label = "TP pied a coulisse"); show(B)
cat("valeurs distinctes pied :", length(unique(tp_pied$y)), " micro :", length(unique(tp_mic$y)), "\n")
## ---- méthode des étendues (AIAG) sur le micromètre
range_method <- function(df) {
  p <- length(unique(df$P)); o <- length(unique(df$O)); r <- nrow(df)/(p*o)
  d2 <- c(`2`=1.128379, `3`=1.692569, `4`=2.058751, `5`=2.325929)
  d3 <- c(`2`=0.852502, `3`=0.888368, `4`=0.879808, `5`=0.864082)
  Rij <- tapply(df$y, list(df$P, df$O), function(v) max(v)-min(v)); Rbb <- mean(Rij)
  sE <- Rbb / d2[as.character(r)]
  xo <- tapply(df$y, df$O, mean); Ro <- max(xo) - min(xo)
  d2s_o <- sqrt(d2[as.character(o)]^2 + d3[as.character(o)]^2)       # une seule étendue (g = 1)
  vO <- max(0, (Ro / d2s_o)^2 - sE^2 / (p * r))
  xp <- tapply(df$y, df$P, mean); Rp <- max(xp) - min(xp)
  ## d2* pour m = p = 10, g = 1 : d2(10) = 3.077505, d3(10) = 0.797478
  sP <- Rp / sqrt(3.077505^2 + 0.797478^2)
  c(Rbb = Rbb, sE = unname(sE), Ro = Ro, d2s_o = unname(d2s_o), sO = sqrt(unname(vO)), Rp = Rp, sP = sP,
    sM = sqrt(unname(sE)^2 + unname(vO)), ndc = sqrt(2) * sP / sqrt(unname(sE)^2 + unname(vO)))
}
print(round(range_method(tp_mic), 5))
cat("K1 = 1/d2(r):", 1/1.128379, 1/1.692569, "  K2 = 1/d2*(o,1):", 1/sqrt(1.128379^2+0.852502^2), 1/sqrt(1.692569^2+0.888368^2),
    "  K3 (10 pièces) :", 1/sqrt(3.077505^2+0.797478^2), "\n")
## ---- données réelles du dépôt
## ---- données réelles du dépôt
m7 <- read.csv2("../mesure_7op_5pieces_3run.csv", dec = "."); names(m7) <- c("P", "O", "y")
show(rr(m7, label = "mesure_7op_5pieces_3run"))
print(round(tapply(m7$y, m7$O, mean), 4)); print(round(tapply(m7$y, m7$O, sd), 4))
x <- rr(subset(m7, O != "op7"), label = "sans op7 (exploratoire)"); show(x)
br <- read.csv2("RR_10pieces_7operateurs_anonyme.csv", dec = ","); names(br) <- c("O", "P", "y")
show(rr(br, label = "RR 10 pieces 7 operateurs 2 essais"))
cat("valeurs extremes :\n"); print(br[br$y >= 39.92, ])
fm <- lmer(y ~ 1 + (1|P) + (1|O), data = br[br$y < 39.92, ]); vc <- as.data.frame(VarCorr(fm))
sP <- sqrt(vc$vcov[vc$grp == "P"]); sM <- sqrt(vc$vcov[vc$grp == "O"] + vc$vcov[vc$grp == "Residual"])
cat("sans valeurs extremes (REML) : sP", sP, " sM", sM, " ndc", sqrt(2) * sP / sM, "\n")
## ---- erreurs de classement (section 6)
LSL <- 11.85; USL <- 12.00; mu <- 11.925
miscl <- function(sP, sM, mu = 11.925) {
  f_in <- function(x) dnorm(x, mu, sP) * (1 - (pnorm(USL, x, sM) - pnorm(LSL, x, sM)))   # X conforme, Y rejeté
  f_out <- function(x) dnorm(x, mu, sP) * (pnorm(USL, x, sM) - pnorm(LSL, x, sM))        # X non conforme, Y accepté
  pin <- pnorm(USL, mu, sP) - pnorm(LSL, mu, sP)
  a <- integrate(f_in, LSL, USL)$value; b <- integrate(f_out, -Inf, LSL)$value + integrate(f_out, USL, Inf)$value
  c(p_conforme = pin, delta = a / pin, beta = b / (1 - pin), rebut_a_tort_ppm = a * 1e6, accepte_a_tort_ppm = b * 1e6)
}
cat("micro :\n"); print(signif(miscl(0.02941, 0.00189), 4))
cat("pied :\n"); print(signif(miscl(0.02984, 0.00975), 4))
cat("Cp vrai/observé micro", 0.15/(6*0.02941), 0.15/(6*0.02947), " pied", 0.15/(6*0.02984), 0.15/(6*0.03139), "\n")
cat("arrondi q/sqrt12 : pied", 0.02/sqrt(12), " micro", 0.001/sqrt(12), "\n")
## cas d'un procédé centré mais capable : sP tel que Cp = 1.33
sPc <- 0.15/(6*1.33); cat("sP pour Cp=1.33 :", sPc, "\n"); print(signif(miscl(sPc, 0.00975), 4)); print(signif(miscl(sPc, 0.00189), 4))
