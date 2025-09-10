##Laboratorio Genética y Bioinformática - Práctica 3##
##Diversidad Genética de Mantas en Isabela y el Continente##

setwd("C:/Users/sebas/Downloads/ej 1")
getwd()
library("hierfstat")
library("ade4")
library("adegenet")

##PARTE I- Diversidad Genética

#Leer matriz de datos para estimar índices de diversidad genética
mantas<-read.csv("Matriz_alélica_mantas_1.csv",header = TRUE,sep = ";")
mantas_diversity<-df2genind(mantas[,-1:-2], ind.names = mantas$X, pop = mantas$Pop,sep = "-", ploidy = 2)

#Estimar la heterocigocidad de cada población
Diversidad_por_locus<-summary(mantas_diversity)
Diversidad_por_locus
Hs_por_poblacion<-Hs(mantas_diversity)
Hs_por_poblacion

#Comparar si los valores de Hs son diferentes entre poblaciones (estadística inferencial)
Hs_ContvsSC<-Hs.test(mantas_diversity[mantas_diversity$pop==1],
                     mantas_diversity[mantas_diversity$pop==2],
                     n.sim = 999,alter = "two-sided")
Hs_ContvsSC


##PARTE II- Distancias Genéticas

#Leer matriz de datos y convertirla a formato para estimar distancias genéticas
mantas2<-read.csv("Matriz_alélica_mantas_2.csv")
mantas2_gen<-as.data.frame(mantas2[,-1])
head(mantas2_gen)

#Cálculo de distancias genéticas

#Distancias de Nei
mantas2_nei<-pairwise.neifst(mantas2_gen,diploid=TRUE)
mantas2_nei
#DIstancias de Weir-Cockram (Fst)
mantas2_FST<-pairwise.WCfst(mantas2_gen, diploid=TRUE)
mantas2_FST

