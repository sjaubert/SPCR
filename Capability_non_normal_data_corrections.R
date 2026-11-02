## Corrections du code de « Analyse de capabilité avec des données non normales »
## S. Jaubert, révision du 2 novembre 2026
## Écrit pour qcc 3.0 (qccGroups, processCapability). Pour qcc 2.x : qcc.groups(x, groupes)
## et process.capability(), mêmes arguments. Version installée : packageVersion("qcc")

library(qcc)
library(MASS)

donnees <- read.csv2("https://sjaubert.github.io/SPCR/device.csv", header = TRUE)
donnees$x <- as.numeric(as.character(donnees$x))

## --- Correction 1 : le regroupement des 100 valeurs en 20 échantillons de 5 ---------
## Avant : rep(1:20, 5) met dans l'échantillon 1 les valeurs 1, 21, 41, 61 et 81.
## Si le fichier est dans l'ordre des prélèvements, un échantillon = 5 valeurs consécutives :
groupes <- rep(1:20, each = 5)
## (garder rep(1:20, 5) seulement si le fichier a été saisi pièce par pièce, colonne par colonne)
d1 <- data.frame(x = donnees$x, groupe = groupes)
diameter <- qccGroups(data = d1, x, groupe)

## --- Correction 2 : ne pas imposer le centre de la carte -------------------------------
## Avant : qcc(diameter, type = "xbar", center = 2)
## La fonction de capabilité calcule Cpk avec le centre de l'objet qcc : avec center = 2,
## il suppose le procédé centré sur la cible et donne Cpk = Cp. La cible se passe
## à processCapability(), pas à qcc().
q1 <- qcc(diameter, type = "xbar")
plot(q1)
mean(donnees$x)          # 1.98757 : le procédé n'est pas centré sur 2

## --- Correction 3 : la capabilité n'a de sens que sur un procédé stable ------------------
## Si la carte montre des points hors limites, le dire avant de calculer les indices :
## les indices décrivent alors un procédé qui n'est pas sous contrôle.
## Ici : avec center = 2, 3 moyennes sortaient sous la LCL ; avec la vraie moyenne,
## il n'en reste qu'une (pour les deux regroupements). Le texte de la page est à reprendre.
processCapability(q1, spec.limits = c(1.95, 2.05), target = 2)

## --- Transformation de Box-Cox (inchangée) ----------------------------------------------
xprime <- donnees$x - 1.8
result <- boxcox(xprime ~ 1, plotit = TRUE, lambda = seq(-4, 4, 0.01))
(lambda <- result$x[which.max(result$y)])        # -2.88
new_donnees <- (xprime^lambda - 1) / lambda
shapiro.test(new_donnees)

New_spec   <- ((c(1.95, 2.05) - 1.8)^lambda - 1) / lambda
New_target <- ((2 - 1.8)^lambda - 1) / lambda

## --- Correction 2 bis : même chose après transformation ---------------------------------
## Avant : qcc(new_diameter, type = "xbar", center = New_center)
## Le centre imposé était la cible transformée (-35.43), alors que la moyenne des
## données transformées vaut -44.65. C'est ce centre imposé, et non le procédé,
## qui produisait Cpk = 0.57.
d2 <- data.frame(x = new_donnees, groupe = groupes)
new_diameter <- qccGroups(data = d2, x, groupe)
new_q <- qcc(new_diameter, type = "xbar")
plot(new_q)
processCapability(new_q, spec.limits = New_spec, target = New_target)

## --- Valeurs attendues (recalculées en Python, sigma = Rbar/d2, d2 = 2.326) ------------
##                         rep(1:20, each = 5)    rep(1:20, 5) (ancien regroupement)
## Avant transformation    sigma 0.0159           sigma 0.0154
##                         Cp 1.05  Cpk 0.79      Cp 1.08  Cpk 0.81
## Après transformation    sigma 10.21            sigma 9.85
##                         Cp 1.03  Cpk 0.86      Cp 1.07  Cpk 0.89
## Pièces hors tolérance observées : 1 sur 100 (2.053).
