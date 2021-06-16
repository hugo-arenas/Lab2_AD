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


#Realizamos analisis de componentes multiple
res.mca <- MCA(tabla, graph = TRUE)

#Para visualizar la proporción de variaciones retenidas por las diferentes dimensiones.
eig.val <- get_eigenvalue(res.mca)

#Para visualizar los porcentajes de inercia explicados por cada dimensión de MCA
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))+ggtitle("Analisis MCA")

#Se muestran las variables con mas contribucion para el estudio

#Dimension 1
fviz_contrib(res.mca, choice = "var", axes = 1, top = 37)+ggtitle("Contribución de variables en Dimension 1")
#Dimension 2
fviz_contrib(res.mca, choice = "var", axes = 2, top = 37)+ggtitle("Contribución de variables en Dimension 2")

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
tabla.scaled <- scale(tabla)

summary(tabla.scaled)

# Se realizan estimaciones por "el método del codo".
print(fviz_nbclust(tabla.scaled, kmeans, method = "wss"))

# Se realizan estimaciones por "el método de la silueta".
print(fviz_nbclust(tabla.scaled, kmeans, method = "silhouette"))

# Se realizan estimaciones por "el método de la brecha estadística".
print(fviz_nbclust(tabla.scaled, kmeans, method = "gap_stat"))

# Tanto el método del codo y el método de la brecha estadística indican que el
# mejor candidato con 3 clusters, sin embargo, el método de la silueta sugiere que
# el mejor candidato con 2 clusters-

# Comprobando el mejor número de clusters por la distancia "euclidean" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters1 <- NbClust(tabla.scaled, distance = "euclidean", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters1))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters1 <- pam(tabla.scaled, k = 2, metric = "euclidean")
print(fviz_cluster(clusters1, data = tabla.scaled, star.plot = TRUE))

dist.eucl = dist(tabla.scaled, method = "euclidean")
print(dist.eucl)
fviz_dist(dist.eucl)

# Comprobando el mejor número de clusters por la distancia "manhattan" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters2 <- NbClust(tabla.scaled, distance = "manhattan", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters2))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters2 <- pam(tabla.scaled, k = 2, metric = "manhattan")
print(fviz_cluster(clusters2, data = tabla.scaled, star.plot = TRUE))

# Comparando ambos clusters, la mejor opción es el obtenido por la distancia de
# "manhattan", dado que este tiene menos datos mezclados entre sus grupos.

#-------------------------------------------------------------------------------

cluster.daisy <- daisy(tabla.scaled, metric = "gower")
cluster.daisy.matrix <- as.matrix(cluster.daisy)
clusters3 <- kmeans(cluster.daisy.matrix, 2)
print(fviz_cluster(clusters3, data = tabla.scaled, star.plot = TRUE))

#-------------------------------------------------------------------------------

dist.manh = dist(tabla.scaled, method = "manhattan")
print(dist.manh)
print(fviz_dist(dist.manh))

gower.dist <- daisy(tabla.scaled, metric = "gower", stand = FALSE)
print(gower.dist)
print(fviz_dist(gower.dist))

gower.matrix <- as.matrix(gower.dist)

sil_width <- c(NA)
for(i in 2:21){  
 pam_fit <- pam(gower.dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}
print(plot(1:21, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width"))
lines(1:21, sil_width)

#Se trabjará con 2 cluster, ya que se sabe que son 2 valores en clases (comestible y venenoso)
#Entonces se quiere ver, cuantos valores pertenecen a cada grupo (cluster)
k<-2


#-------------------------------------------------------------------------------

tabla$clusters <- clusters2$clustering
print(summary(tabla[tabla["clusters"] == 1, ]))
print(summary(tabla[tabla["clusters"] == 2, ]))
print(summary(tabla["clusters"]))

# Se obtienen los medioides reales respecto a cada variable por cluster.
med.class <- clusters2$medoids[,1] * attr(tabla.scaled, 'scaled:scale')[1] + attr(tabla.scaled, 'scaled:center')[1]
med.age <- clusters2$medoids[,2] * attr(tabla.scaled, 'scaled:scale')[2] + attr(tabla.scaled, 'scaled:center')[2]
med.menopause <- clusters2$medoids[,3] * attr(tabla.scaled, 'scaled:scale')[3] + attr(tabla.scaled, 'scaled:center')[3]
med.tumor.size <- clusters2$medoids[,4] * attr(tabla.scaled, 'scaled:scale')[4] + attr(tabla.scaled, 'scaled:center')[4]
med.inv.nodes <- clusters2$medoids[,5] * attr(tabla.scaled, 'scaled:scale')[5] + attr(tabla.scaled, 'scaled:center')[5]
med.node.caps <- clusters2$medoids[,6] * attr(tabla.scaled, 'scaled:scale')[6] + attr(tabla.scaled, 'scaled:center')[6]
med.deg.malig <- clusters2$medoids[,7] * attr(tabla.scaled, 'scaled:scale')[7] + attr(tabla.scaled, 'scaled:center')[7]
med.breast <- clusters2$medoids[,8] * attr(tabla.scaled, 'scaled:scale')[8] + attr(tabla.scaled, 'scaled:center')[8]
med.breast.quad <- clusters2$medoids[,9] * attr(tabla.scaled, 'scaled:scale')[9] + attr(tabla.scaled, 'scaled:center')[9]
med.irradiat <- clusters2$medoids[,10] * attr(tabla.scaled, 'scaled:scale')[10] + attr(tabla.scaled, 'scaled:center')[10]

# Se crea una tabla que muestra los medioides de cada variables según los grupos
# de cluster
med.tabla <- data.frame(med.class, med.age, med.menopause, med.tumor.size,
                        med.inv.nodes, med.node.caps, med.deg.malig, med.breast,
                        med.breast.quad, med.irradiat)



