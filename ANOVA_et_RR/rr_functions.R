## Fonctions pour l'étude R&R (modèle aléatoire croisé, plan équilibré)
## Références : Burdick, Borror, Montgomery (2005), chap. 3 et 5.

mls_constants <- function(p, o, r, alpha = 0.05) {
  a2 <- alpha / 2
  c(G1 = 1 - qf(a2, Inf, p - 1), G2 = 1 - qf(a2, Inf, o - 1),
    G3 = 1 - qf(a2, Inf, (p - 1) * (o - 1)), G4 = 1 - qf(a2, Inf, p * o * (r - 1)),
    H1 = qf(1 - a2, Inf, p - 1) - 1, H2 = qf(1 - a2, Inf, o - 1) - 1,
    H3 = qf(1 - a2, Inf, (p - 1) * (o - 1)) - 1, H4 = qf(1 - a2, Inf, p * o * (r - 1)) - 1,
    F1 = qf(1 - a2, p - 1, (p - 1) * (o - 1)), F2 = qf(a2, p - 1, (p - 1) * (o - 1)),
    F3 = qf(1 - a2, p - 1, o - 1), F4 = qf(a2, p - 1, o - 1))
}

## Intervalles MLS, modèle avec interaction (Burdick et al., éq. 3.3 à 3.5)
mls_interaction <- function(SP, SO, SPO, SE, p, o, r, alpha = 0.05) {
  k <- as.list(mls_constants(p, o, r, alpha))
  with(k, {
    G13 <- ((F1 - 1)^2 - G1^2 * F1^2 - H3^2) / F1
    H13 <- ((1 - F2)^2 - H1^2 * F2^2 - G3^2) / F2
    gP <- (SP - SPO) / (o * r)
    gM <- (SO + (p - 1) * SPO + p * (r - 1) * SE) / (p * r)
    VLP <- G1^2 * SP^2 + H3^2 * SPO^2 + G13 * SP * SPO
    VUP <- H1^2 * SP^2 + G3^2 * SPO^2 + H13 * SP * SPO
    VLM <- G2^2 * SO^2 + G3^2 * (p - 1)^2 * SPO^2 + G4^2 * p^2 * (r - 1)^2 * SE^2
    VUM <- H2^2 * SO^2 + H3^2 * (p - 1)^2 * SPO^2 + H4^2 * p^2 * (r - 1)^2 * SE^2
    LR <- p * (1 - G1) * (SP - F1 * SPO) / (p * o * (r - 1) * SE + o * (1 - G1) * F3 * SO + o * (p - 1) * SPO)
    UR <- p * (1 + H1) * (SP - F2 * SPO) / (p * o * (r - 1) * SE + o * (1 + H1) * F4 * SO + o * (p - 1) * SPO)
    data.frame(parametre = c("gamma_P", "gamma_M", "gamma_R"),
               estimation = c(gP, gM, gP / gM),
               L = pmax(0, c(gP - sqrt(VLP) / (o * r), gM - sqrt(VLM) / (p * r), LR)),
               U = c(gP + sqrt(VUP) / (o * r), gM + sqrt(VUM) / (p * r), UR),
               G13 = G13, H13 = H13)
  })
}

## Intervalles généralisés (GPQ, Burdick et al., tableau 3.8)
gci_interaction <- function(SP, SO, SPO, SE, p, o, r, alpha = 0.05, N = 1e5, seed = 1) {
  set.seed(seed)
  W1 <- rchisq(N, p - 1); W2 <- rchisq(N, o - 1); W3 <- rchisq(N, (p - 1) * (o - 1)); W4 <- rchisq(N, p * o * (r - 1))
  gP <- pmax(0, (p - 1) * SP / (o * r * W1) - (p - 1) * (o - 1) * SPO / (o * r * W3))
  gM <- (o - 1) * SO / (p * r * W2) + (p - 1)^2 * (o - 1) * SPO / (p * r * W3) + p * o * (r - 1)^2 * SE / (r * W4)
  q <- function(x) quantile(x, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  rbind(gamma_P = q(gP), gamma_M = q(gM), gamma_R = q(gP / gM))
}

## Intervalles MLS, modèle sans interaction (Burdick et al., éq. 5.3 à 5.5)
mls_additive <- function(SP, SO, SE, p, o, r, alpha = 0.05) {
  a2 <- alpha / 2; nE <- p * o * r - p - o + 1
  G1 <- 1 - qf(a2, Inf, p - 1); G2 <- 1 - qf(a2, Inf, o - 1); G3 <- 1 - qf(a2, Inf, nE)
  H1 <- qf(1 - a2, Inf, p - 1) - 1; H2 <- qf(1 - a2, Inf, o - 1) - 1; H3 <- qf(1 - a2, Inf, nE) - 1
  F1 <- qf(1 - a2, p - 1, nE); F2 <- qf(a2, p - 1, nE); F3 <- qf(1 - a2, p - 1, o - 1); F4 <- qf(a2, p - 1, o - 1)
  G13 <- ((F1 - 1)^2 - G1^2 * F1^2 - H3^2) / F1; H13 <- ((1 - F2)^2 - H1^2 * F2^2 - G3^2) / F2
  gP <- (SP - SE) / (o * r); gM <- (SO + (p * r - 1) * SE) / (p * r)
  VLP <- G1^2 * SP^2 + H3^2 * SE^2 + G13 * SP * SE; VUP <- H1^2 * SP^2 + G3^2 * SE^2 + H13 * SP * SE
  VLM <- G2^2 * SO^2 + G3^2 * (p * r - 1)^2 * SE^2; VUM <- H2^2 * SO^2 + H3^2 * (p * r - 1)^2 * SE^2
  LR <- p * ((1 - G1) * SP^2 - SP * SE + (F1 - (1 - G1) * F1^2) * SE^2) / (o * ((p * r - 1) * SP * SE + (1 - G1) * F3 * SP * SO))
  UR <- p * ((1 + H1) * SP^2 - SP * SE + (F2 - (1 + H1) * F2^2) * SE^2) / (o * ((p * r - 1) * SP * SE + (1 + H1) * F4 * SP * SO))
  data.frame(parametre = c("gamma_P", "gamma_M", "gamma_R"), estimation = c(gP, gM, gP / gM),
             L = pmax(0, c(gP - sqrt(VLP) / (o * r), gM - sqrt(VLM) / (p * r), LR)),
             U = c(gP + sqrt(VUP) / (o * r), gM + sqrt(VUM) / (p * r), UR))
}
