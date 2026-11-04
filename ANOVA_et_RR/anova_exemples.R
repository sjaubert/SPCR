## anova_exemples.R : calculs du document « L'analyse de la variance à plusieurs facteurs »
## S. Jaubert, novembre 2026. R >= 4.0, paquet car pour la section 7.
## À lancer depuis le dossier ANOVA_et_RR du dépôt SPCR (rugo_vib.csv est dans le dossier parent).
options(digits = 7)
## ---------- Exemple 1 : batteries (2 voltmètres x 3 batteries x 3 essais)
voltage <- c(1.4727,1.4206,1.4754,1.5083,1.5739,1.4341,1.5517,1.5483,1.4614,1.3337,1.6078,1.4767,1.4066,1.5951,1.8419,1.7087,1.8259,1.5444)
bat <- data.frame(A = factor(rep(1:2, each = 9)), B = factor(rep(rep(1:3, each = 3), 2)), y = voltage)
a <- 2; b <- 3; r <- 3
m <- mean(bat$y); mi <- tapply(bat$y, bat$A, mean); mj <- tapply(bat$y, bat$B, mean); mij <- tapply(bat$y, list(bat$A, bat$B), mean)
SCT <- sum((bat$y - m)^2); SCA <- b*r*sum((mi-m)^2); SCB <- a*r*sum((mj-m)^2)
SCAB <- r*sum((mij - outer(mi, rep(1,b)) - outer(rep(1,a), mj) + m)^2)
SCR <- sum((bat$y - mij[cbind(bat$A, bat$B)])^2)
cat("Ex1 moyennes:\n"); print(round(mi,5)); print(round(mj,5)); print(round(mij,5)); cat("m=",m,"\n")
cat(sprintf("SCT %.6f SCA %.6f SCB %.6f SCAB %.6f SCR %.6f somme %.6f\n", SCT, SCA, SCB, SCAB, SCR, SCA+SCB+SCAB+SCR))
fit <- aov(y ~ A*B, data = bat); print(summary(fit))
MS <- c(A = SCA/(a-1), B = SCB/(b-1), AB = SCAB/((a-1)*(b-1)), R = SCR/(a*b*(r-1)))
cat("MS:", MS, "\n")
cat("Modèle fixe F :", MS[1:3]/MS[4], " p :", pf(MS[1:3]/MS[4], c(a-1,b-1,(a-1)*(b-1)), a*b*(r-1), lower.tail=FALSE), "\n")
cat("Modèle aléatoire F_A=MS_A/MS_AB", MS["A"]/MS["AB"], "p", pf(MS["A"]/MS["AB"], a-1, (a-1)*(b-1), lower.tail=FALSE),
    " F_B", MS["B"]/MS["AB"], "p", pf(MS["B"]/MS["AB"], b-1, (a-1)*(b-1), lower.tail=FALSE), "\n")
cat("Seuils F 5% :", qf(0.95, 1, 12), qf(0.95, 2, 12), qf(0.95,1,2), qf(0.95,2,2), "\n")
## ---------- Exemple 2 : flash mesure, 7 étudiants x 5 butées, sans répétition
X <- matrix(c(11.92,11.94,11.98,11.96,11.96, 11.94,11.92,11.92,11.94,11.94, 11.92,11.92,11.94,11.94,11.96,
              11.94,11.90,11.92,11.96,11.96, 11.94,11.90,11.94,11.92,11.90, 11.94,11.92,11.94,11.96,11.94,
              11.94,11.92,11.96,12.00,11.91), nrow = 7, byrow = TRUE)
fl <- data.frame(O = factor(rep(1:7, 5)), P = factor(rep(1:5, each = 7)), y = as.vector(X))
print(summary(aov(y ~ O + P, data = fl)))
## test de Tukey à 1 ddl
res_tukey <- function(X) {
  a <- nrow(X); b <- ncol(X); m <- mean(X); ai <- rowMeans(X) - m; bj <- colMeans(X) - m
  SCN <- (sum(X * outer(ai, bj)))^2 / (sum(ai^2) * sum(bj^2))
  SCres <- sum((X - outer(rowMeans(X), rep(1,b)) - outer(rep(1,a), colMeans(X)) + m)^2)
  dfr <- (a-1)*(b-1) - 1; Fn <- SCN / ((SCres - SCN)/dfr)
  c(SCN = SCN, SCres = SCres, F = Fn, df2 = dfr, p = pf(Fn, 1, dfr, lower.tail = FALSE))
}
print(res_tukey(X))
## vérification par régression (méthode équivalente)
f0 <- lm(y ~ O + P, data = fl); fl$q <- fitted(f0)^2; print(anova(lm(y ~ O + P + q, data = fl)))
X2 <- X; X2[7,5] <- 11.92  # hypothèse : 11,91 était 11,92
print(summary(aov(y ~ O + P, data = transform(fl, y = as.vector(X2)))))

## ---------- Exemple 3 : plan d'usinage 3^3 sans répétition
d <- read.csv2("../rugo_vib.csv"); d$D<-factor(d$D); d$V<-factor(d$V); d$f<-factor(d$f)
cat("moyenne Ra", mean(d$Ra), "\n")
full <- aov(Ra ~ D*V*f, data=d); s <- summary(full)[[1]]; print(s)
cat("SCT", sum((d$Ra-mean(d$Ra))^2), " somme SC", sum(s[,"Sum Sq"]), "\n")
m2 <- aov(Ra ~ (D+V+f)^2, data=d); print(summary(m2))
## regroupement des interactions non significatives à 5 % (D:V et V:f) avec le résidu
m3 <- aov(Ra ~ D+V+f+D:f, data=d); print(summary(m3))
print(model.tables(aov(Ra ~ D+V+f+D:f, data=d), "means")$tables$`D:f`)
print(round(tapply(d$Ra, list(d$D, d$f), mean),3))
## résidus
r2 <- residuals(m2); print(shapiro.test(r2))
## La simulation sur l'interaction ABC est dans sim_abc.R
## ---- trois facteurs aléatoires (Burdick tableau 8.18) : quasi-F pour Parts
S <- c(P=102.8, O=13.08, T=29.23, PO=1.074, PT=1.884, OT=2.932, POT=0.4183, E=0.1394)
n <- c(P=8, O=2, T=2, PO=16, PT=16, OT=4, POT=32, E=81)
num <- S["P"] + S["POT"]; den <- S["PO"] + S["PT"]
nu1 <- num^2/(S["P"]^2/n["P"] + S["POT"]^2/n["POT"]); nu2 <- den^2/(S["PO"]^2/n["PO"] + S["PT"]^2/n["PT"])
Fq <- num/den; cat("quasi-F parts", Fq, "ddl", nu1, nu2, "p", pf(Fq, nu1, nu2, lower.tail=FALSE), "\n")
## variante Fq = MS_P / (MS_PO + MS_PT - MS_POT)
den2 <- S["PO"] + S["PT"] - S["POT"]; nu2b <- den2^2/(S["PO"]^2/n["PO"] + S["PT"]^2/n["PT"] + S["POT"]^2/n["POT"])
cat("quasi-F' ", S["P"]/den2, "ddl", n["P"], nu2b, "\n")

## ---------- Section 7 : plan déséquilibré (types I, II, III)
d <- bat[-c(1, 2), ]; print(table(d$A, d$B))
print(anova(lm(y ~ A*B, d))); print(anova(lm(y ~ B*A, d)))
print(car::Anova(lm(y ~ A*B, d), type = 2))
op <- options(contrasts = c("contr.sum", "contr.poly")); print(car::Anova(lm(y ~ A*B, d), type = 3)); options(op)
