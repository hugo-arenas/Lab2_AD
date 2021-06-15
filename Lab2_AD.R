#Laboratorio 2 - Análisis de Datos
#Integrantes:
#            -Hugo Arenas
#            -Juan Arredondo

library(cluster)
library(fpc)
library(NbClust)
library(factoextra)
library(FactoMineR)
library(Rtsne) 

library(ez)
library(dplyr)
library(ggpubr)
library(knitr)
library(tidyr)
library(car)
library(lmtest)
library(vcd)

#Los nombres originales de las variables, una breve explicación y los tipos de los datos:

#age (rango de edad): 10-19, 20-29, 30-39, 40-49,.
#menopause (momento de la menopausia): lt40, ge40, premeno.
#tumor-size (tamaño del tumor extirpado en mm): 0-4, 5-9, 10-14, .
#inv-nodes (una métrica de presencia de células cancerosas en los nodos linfáticos): 0-2, 3-5, 6-8, 9-11,.
#node-caps (evidencia de que células cancerosas atravesaron la cápsula de los nódulos linfáticos): yes, no
#deg-malig (grado histológico del tumor: bajo, intermedio, alto): 1, 2, 3.
#breast (mama afectada): left, right.
#breast-quad (cuadrante de la mama): left-up, left-low, right-up, right-low, central.
#irradiat (radioterapia): yes, no.
#Class (clase) Indica recurrencia, es la variable a predecir (no-recurrencia: 201 casos, recurrencia: 85 casos)


# Leemos los datos
dirstudio <- dirname(rstudioapi::getSourceEditorContext()$path)
filename <- "breast-cancer.data"
file <- file.path(dirstudio, filename)

#Se definen los nombres de las columnas, estos son los mismos de provistos por la base de datos
columns <- c("class", 
             "age", 
             "menopause", 
             "tumor.size", 
             "inv.nodes", 
             "node.caps",
             "deg.malig",
             "breast",
             "breast.quad",
             "irradiat")

tabla <- read.csv(file, col.names = columns)
tabla$class <- unclass(as.factor(tabla$class))
tabla$age <- unclass(as.factor(tabla$age)) 
tabla$menopause <- unclass(as.factor(tabla$menopause)) 
tabla$tumor.size <- unclass(as.factor(tabla$tumor.size))
tabla$inv.nodes <- unclass(as.factor(tabla$inv.nodes)) 
tabla$node.caps <- unclass(as.factor(tabla$node.caps))
tabla$deg.malig <- unclass(as.factor(tabla$deg.malig))
tabla$breast <- unclass(as.factor(tabla$breast)) 
tabla$breast.quad <- unclass(as.factor(tabla$breast.quad)) 
tabla$irradiat <- unclass(as.factor(tabla$irradiat))

#Lo primero es ver el rango de las variables, algunas medidas
#de tendencia central para las variables continuas, recuentos 
#para las categóricas y presencia de valores faltantes.
summary(tabla)

#Se sacan los datos nulos del datagrama
bool.values <- tabla$node.caps=='1'
tabla <- tabla[!bool.values,]

bool.values <- tabla$breast.quad =='1'
tabla <- tabla[!bool.values,]

#Se pone a escala.
tabla <- scale(tabla)

summary(tabla)

# Se realizan estimaciones por "el método del codo".
print(fviz_nbclust(tabla, kmeans, method = "wss"))

# Se realizan estimaciones por "el método de la silueta".
print(fviz_nbclust(tabla, kmeans, method = "silhouette"))

# Se realizan estimaciones por "el método de la brecha estadística".
print(fviz_nbclust(tabla, kmeans, method = "gap_stat"))

# Tanto el método del codo y el método de la brecha estadística indican que el
# mejor candidato con 3 clusters, sin embargo, el método de la silueta sugiere que
# el mejor candidato con 2 clusters-

# Comprobando el mejor número de clusters por la distancia "euclidean" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters1 <- NbClust(tabla, distance = "euclidean", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters1))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters1 <- pam(tabla, k = 2, metric = "euclidean")
print(fviz_cluster(clusters1, data = tabla, star.plot = TRUE))

# Comprobando el mejor número de clusters por la distancia "manhattan" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters2 <- NbClust(tabla, distance = "manhattan", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters2))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters2 <- pam(tabla, k = 2, metric = "manhattan")
print(fviz_cluster(clusters2, data = tabla, star.plot = TRUE))

# Comparando ambos clusters, la mejor opción es el obtenido por la distancia de
# "manhattan", dado que este tiene menos datos mezclados entre sus grupos.


gower.dist <- daisy(tabla, metric = c("gower"))
print(fviz_nbclust(gower.dist))
gower.matrix = as.matrix(gower.dist)
tabla[
  which(gower.matrix == min(gower.matrix[gower.matrix != min(gower.matrix)]),
        arr.ind = TRUE)[1, ], ]

# Calculate silhouette width for many k using PAM

sil_width <- c(NA)

for(i in 2:10){
  
  pam_fit <- pam(gower.dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot sihouette width (higher is better)

plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

pam_fit <- pam(gower.dist, diss = TRUE, k = 2)

tsne_obj <- Rtsne(gower.dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = tabla$class)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))
#clusters3 <- pam(gower.dist, k = 2, diss = TRUE)
#print(fviz_cluster(clusters3, data = tabla, star.plot = TRUE))
