## TP n°2 ACF ## 
## PAOLO _ LEROY _ GIBOUDEAU _ BATTISTINI ##

library(openxlsx)
library(FactoMineR)
library(ggplot2)
library(grid)
source("ACP_param.r")
source("quality_ACP.r")


#### Partie 2 : Base de données "TP_AFC_majeur1718_travail.xlsx" ####

# Question 1

travail <- read.xlsx("TP_AFC_majeur1718_travail.xlsx")

travail_propres <- travail[!(travail$temps.travail == 99.0 | travail$Fonction == 0),]

fonction<-travail_propres$Fonction
temps_travail<-travail_propres$temps.travail

fonction<-na.omit(fonction)
temps_travail<-na.omit(temps_travail)

fonction_unique<-unique(fonction)
temps_travail_unique<-unique(temps_travail)

fonction_unique<-na.omit(fonction_unique)
temps_travail_unique<-na.omit(temps_travail_unique)

fonction_unique<-fonction_unique[fonction_unique!=" "]

fonction_unique<-sort(fonction_unique)
temps_travail_unique<-sort(temps_travail_unique)

#Table de contingence

table_contingence <- matrix(0, nrow = length(fonction_unique), ncol = length(temps_travail_unique))

for (i in fonction_unique){
  for (j in temps_travail_unique){
    count<-0
    for (k in 1:length(fonction)){
      if ((fonction[k]==i) & (temps_travail[k]==j)){
        count<-count+1
        
      }}
    table_contingence[which(fonction_unique==i), which(temps_travail_unique==j)]<-count
  }
}

rownames(table_contingence) <- c("1", "2", "3", "4", "5", "6", "7")
colnames(table_contingence) <- c("2", "4", "8", "16")

# Calculer le nombre total d'observations
nb_obs <- sum(table_contingence)

# Calculer le tableau des fréquences relatives
frequence <- table_contingence / nb_obs

# Conversion du tableau de fréquence en matrice

frequence_mat <- as.matrix(frequence)

# Ajout des marges au tableau de fréquences relatives
frequence_marge <- as.data.frame(addmargins(frequence_mat))

# Calcul du profil ligne moyen
D_n <- diag(margin.table(frequence_mat,1))
Lignes <- solve(D_n) %*% frequence_mat


# Calcul des profils colonnes (matriciel)
D_p <- diag(margin.table(frequence_mat,2))
colonnes <- solve(D_p) %*% t(frequence_mat)

#Question 2

# Matrice dont on cherche les vecteurs propres dans l'analyse directe
S2 <- t(frequence_mat) %*% solve(D_n) %*% frequence_mat %*% solve(D_p)

valeur_propres <- eigen(S2)$values
vecteurs_propres <- eigen(S2)$vectors

comp <- Lignes %*% vecteurs_propres
Comp <- comp

#Question 3 et 4

#Affichage de la projection 

plot.new()
plot(Comp[,2], Comp[,1],  xlab="C1", ylab="C2", ylim=c(-0.1, 0.8), xlim = c(-0.55, 0.4), col='Blue')
text(Comp[,2], Comp[,1]+ 0.03, labels = c("1", "2", "3", "4", "5", "6", "7"), col='Blue')
abline(0, 0, col='Red')
abline(0, 1000, col='Red')

#Analyse des valeurs propres en cascade

quality_ACP(valeur_propres, Comp)


# Question 5

summary(table_contingence)

res.test.chi2<- chisq.test(table_contingence)
round(res.test.chi2$expected,1)
round(res.test.chi2$residuals^2, 2)
round(100 * res.test.chi2$residuals^2 / res.test.chi2$stat, 2)


# Projection du profil colonne 

# Matrice dont on cherche les vecteurs propres dans l'analyse directe
TT <- frequence_mat %*% solve(D_p) %*% t(frequence_mat) %*% solve(D_n)

valeur_propres_2 <- eigen(TT)$values
vecteurs_propres_2 <- eigen(TT)$vectors

comp_2 <- colonnes %*% vecteurs_propres_2
Comp_2 <- comp_2

#Analyse des valeurs propres en cascade

quality_ACP(valeur_propres_2, Comp_2)


# Affichage des projections du profil des lignes et du profil des colonnes 
plot.new()
plot(Comp_2[,1], Comp_2[,2], col='Green', ylim=c(-0.3, 0.3), xlim = c(-0.6, 0.05))
abline(0, 0, col='Red')
abline(0, 1000, col='Red')
text(Comp_2[,1], (Comp_2[,2])+0.03, labels = c("2", "4", "8", "16"), col='Green')
grid()

### Partie avec FactoMineR

# Effectuer l'Analyse Factorielle des Correspondances (AFC)
afc_result_pb <- CA(table_contingence)

#### Partie 1 : Exemple de cours ####

### Partie codée main

# Créer la matrice de données à partir du tableau de contingence
kij <- matrix(c(10, 0, 0, 0, 8, 2, 0, 4, 6), nrow = 3, byrow = TRUE)
rownames(kij) <- c("Sucre", "Acide", "Amer")
colnames(kij) <- c("Perçu sucré", "Perçu acide", "Perçu amer")

# Calculer le nombre total d'observations
kpp <- sum(kij)

# Calculer le tableau des fréquences relatives
fij <- kij / kpp

# Conversion du tableau de fréquence en matrice

fij_mat <- as.matrix(fij)

# Ajout des marges au tableau de fréquences relatives
fij_marge <- as.data.frame(addmargins(fij_mat))

# Calcul du profil ligne moyen
Dn <- diag(margin.table(fij_mat,1))
L <- solve(Dn) %*% fij_mat


# Calcul des profils colonnes (matriciel)
Dp <- diag(margin.table(fij_mat,2))
C <- t(solve(Dp) %*% t(fij_mat))
rownames(C) <- colnames(fij_mat)


# Matrice dont on cherche les vecteurs propres dans l'analyse directe
S <- t(fij_mat) %*% solve(Dn) %*% fij_mat %*% solve(Dp)

Val_p <- eigen(S)$values
Vec_p <- eigen(S)$vectors

comp <- L %*% Vec_p
Comp <- comp[,2:3]

#Analyse des valeurs propres en cascade

quality_ACP(Val_p[2:3], Comp)

# Matrice dont on cherche les vecteurs propres dans l'analyse indirecte
TT <- fij_mat %*% solve(Dp) %*% t(fij_mat) %*% solve(Dn)

Val_p_2 <- eigen(TT)$values
Vec_p_2 <- eigen(TT)$vectors

comp_2 <- C %*% Vec_p_2
Comp_2 <- comp_2[,2:3]

#Analyse des valeurs propres en cascade

quality_ACP(Val_p_2[2:3], Comp_2)

point_names_l <- c("Sucre", "Acide", "Amer")
point_names_c <- c("Perçu sucré", "Perçu acide", "Perçu amer")


# Affichage des projections du profil des lignes et du profil des colonnes 
plot.new()
plot(Comp[,1], Comp[,2],  xlab="C1", ylab="C2", ylim=c(-0.45, 0.45), xlim = c(-0.07,1.1), col='Blue')
text(Comp[,1], Comp[,2]+ 0.03, labels = c("Sucre", "Acide", "Amer"), col='Black')
points(Comp_2[,1], Comp_2[,2], col='Green', xlab="C1", ylab="C2", ylim=c(-0.45, 0.45), xlim = c(-0.07,1.1))
abline(0, 0, col='Red')
abline(0, 1000, col='Red')
text(Comp_2[,1], Comp_2[,2]+ 0.06, labels = c("Perçu sucré", "Perçu acide", "Perçu amer"), col='Black', ps=30)
grid()

### Partie avec FactoMineR

# Effectuer l'Analyse Factorielle des Correspondances (AFC)
afc_result <- CA(kij)

