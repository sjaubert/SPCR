## student_khi2.R : calculs du document « Pourquoi Student, pourquoi le khi-deux »
## S. Jaubert, novembre 2026. R >= 4.0, sans paquet.
options(digits = 6)
## ---- 1. exemple : premier prélèvement de 5 butées (diametre_butee.txt, dépôt SPCR)
x <- c(11.87, 11.86, 11.84, 11.88, 11.87); n <- length(x); nu <- n - 1
m <- mean(x); s <- sd(x); cat("n", n, "moyenne", m, "s", s, "\n")
tq <- qt(0.975, nu); zq <- qnorm(0.975)
cat("t", tq, " IC moyenne", m - tq*s/sqrt(n), m + tq*s/sqrt(n), " demi-largeur", tq*s/sqrt(n), "\n")
cat("faux IC z avec s", m - zq*s/sqrt(n), m + zq*s/sqrt(n), " demi-largeur", zq*s/sqrt(n), "\n")
q <- qchisq(c(0.025, 0.975), nu); cat("chi2 quantiles", q, "\n")
cat("IC variance", nu*s^2/q[2], nu*s^2/q[1], " IC sigma", sqrt(nu*s^2/q[2]), sqrt(nu*s^2/q[1]), "\n")
cat("multiplicateurs sigma", sqrt(nu/q[2]), sqrt(nu/q[1]), "\n")
cat("P(|T|<=1.96) nu=4 :", 2*pt(zq, 4) - 1, "\n")
## les 10 prélèvements : moyennes
d <- matrix(c(11.87,11.86,11.84,11.88,11.87, 11.87,11.85,11.86,11.84,11.89, 11.91,11.9,11.86,11.88,11.87,
 11.89,11.91,11.89,11.88,11.92, 11.91,11.92,11.93,11.89,11.9, 11.85,11.91,11.92,11.93,11.9,
 11.91,11.92,11.95,11.94,11.96, 11.96,11.94,11.95,11.97,11.94, 11.97,11.96,11.96,11.95,11.94,
 12.01,12.04,12.03,11.96,11.99), 10, byrow = TRUE)
cat("moyennes des prélèvements", round(rowMeans(d), 3), "\n")
cat("s des prélèvements", round(apply(d, 1, sd), 4), " s pooled", sqrt(mean(apply(d,1,var))), " s global", sd(as.vector(d)), "\n")
summary(aov(v ~ g, data.frame(v = as.vector(t(d)), g = factor(rep(1:10, each = 5))))) |> print()
## ---- 2. quantiles
nus <- c(1, 2, 4, 9, 19, 29, 49, 99)
print(round(rbind(nu = nus, t975 = qt(0.975, nus), chi025 = qchisq(0.025, nus), chi975 = qchisq(0.975, nus)), 3))
cat("z975", zq, "\n")
cat("variance de T nu = 4, 9, 29 :", 4/2, 9/7, 29/27, "\n")
cat("P(T>3) nu=4 :", pt(3, 4, lower.tail = FALSE), " normale :", pnorm(3, lower.tail = FALSE), "\n")
## ---- 3. simulation : 40 intervalles (figure)
set.seed(2026); K <- 40; nn <- 10
ic <- t(replicate(K, { y <- rnorm(nn, 10, 1); mean(y) + c(-1, 1) * qt(0.975, nn-1) * sd(y)/sqrt(nn) }))
cover <- ic[,1] <= 10 & ic[,2] >= 10; cat("intervalles qui contiennent mu :", sum(cover), "/", K, "\n")
write.table(round(cbind(i = 1:K, L = ic[,1], U = ic[,2], ok = as.integer(cover)), 3), "ic40.dat", row.names = FALSE, quote = FALSE)
## ---- 4. robustesse (taux de couverture à 95 %)
set.seed(1); N <- 20000
lois <- list(normale = function(n) rnorm(n), uniforme = function(n) runif(n), exponentielle = function(n) rexp(n), student5 = function(n) rt(n, 5))
mu <- c(0, 0.5, 1, 0); sig <- c(1, sqrt(1/12), 1, sqrt(5/3)); kurt <- c(0, -1.2, 6, 6)
res <- NULL
for (L in seq_along(lois)) for (n in c(5, 10, 30, 100)) {
  cm <- cs <- cz <- 0
  for (k in 1:N) { y <- lois[[L]](n); mm <- mean(y); ss <- sd(y)
    h <- qt(0.975, n-1)*ss/sqrt(n); cm <- cm + (abs(mm - mu[L]) <= h)
    cz <- cz + (abs(mm - mu[L]) <= qnorm(0.975)*ss/sqrt(n))
    v <- (n-1)*ss^2; cs <- cs + (v/qchisq(0.975, n-1) <= sig[L]^2 & sig[L]^2 <= v/qchisq(0.025, n-1)) }
  res <- rbind(res, data.frame(loi = names(lois)[L], n = n, student = cm/N, z_avec_s = cz/N, khi2 = cs/N))
}
print(res)
cat("couverture asymptotique khi2 si kurtosis kappa :", round(2*pnorm(qnorm(0.975)/sqrt(1 + kurt/2)) - 1, 3), "\n")
## ---- 5. taille d'échantillon
nmean <- function(e) { n <- 2; while (qt(0.975, n-1)/sqrt(n) > e) n <- n + 1; n }   # e = demi-largeur / sigma
cat("n pour E/sigma = 1, 0.5, 0.25, 0.1 :", sapply(c(1, 0.5, 0.25, 0.1), nmean), "\n")
ns <- c(5, 10, 20, 30, 50, 100, 200, 500)
print(round(rbind(n = ns, bas = sqrt((ns-1)/qchisq(0.975, ns-1)), haut = sqrt((ns-1)/qchisq(0.025, ns-1))), 3))
## ---- 6. capabilité
for (n in c(25, 50, 100, 125)) cat("n", n, " Cp=1.33 -> IC", round(1.33*sqrt(qchisq(c(0.025, 0.975), n-1)/(n-1)), 3),
  " borne inf. unilatérale 95 %", round(1.33*sqrt(qchisq(0.05, n-1)/(n-1)), 3), "\n")
biais <- function(n) { nu <- n-1; sqrt(nu)*gamma((nu-1)/2)/(sqrt(2)*gamma(nu/2)) }
cat("E(Cp chapeau)/Cp n = 5, 10, 25, 50, 125 :", round(sapply(c(5,10,25,50,125), biais), 4), "\n")
c4 <- function(n) sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2); cat("c4 :", round(sapply(c(5,10,25,50,125), c4), 4), "\n")
set.seed(3); cat("contrôle simulation n=10 :", mean(replicate(2e5, 1/sd(rnorm(10)))), "\n")
## exemple Cp prélèvement 1
cat("Cp prélèvement 1 :", 0.15/(6*s), " IC", 0.15/(6*s)*sqrt(q/nu), "\n")
## ---- 7. loi F
s1 <- sd(d[1,]); s2 <- sd(d[10,]); Fobs <- s2^2/s1^2
cat("F", Fobs, " p bilatérale", 2*min(pf(Fobs, 4, 4), pf(Fobs, 4, 4, lower.tail = FALSE)), " quantiles", qf(c(0.025, 0.975), 4, 4), "\n")
cat("IC rapport des variances", Fobs/qf(0.975, 4, 4), Fobs/qf(0.025, 4, 4), "\n")
cat("t^2 = F(1,nu) : qt(0.975,4)^2", qt(0.975,4)^2, " qf(0.95,1,4)", qf(0.95,1,4), "\n")
print(var.test(d[10,], d[1,]))
print(t.test(x)$conf.int)
