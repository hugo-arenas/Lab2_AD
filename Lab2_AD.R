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
print(fviz_nbclust(tabla.scaled, kmeans, method = "wss")+
        labs(title = "Método del codo"))

# Se realizan estimaciones por "el método de la silueta".
print(fviz_nbclust(tabla.scaled, kmeans, method = "silhouette")+
        labs(title = "Método de la Silueta"))

# Se realizan estimaciones por "el método de la brecha estadística".
print(fviz_nbclust(tabla.scaled, kmeans, method = "gap_stat")+
        labs(title = "Método de la brecha estadística"))

# Se realizan estimaciones por "el método de la brecha estadística".
print(fviz_nbclust(tabla.scaled, kmeans, method = "gap_stat")+
        labs(title = "Método de la brecha estadística"))

# Tanto el método del codo y el método de la brecha estadística indican que el
# mejor candidato con 3 clusters, sin embargo, el método de la silueta sugiere que
# el mejor candidato con 2 clusters-

# Comprobando el mejor número de clusters por la distancia "euclidean" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters1 <- NbClust(tabla.scaled, distance = "euclidean", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters1) +
        labs(title = "K óptimo con distancia Euclediana"))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters1 <- pam(tabla.scaled, k = 2, metric = "euclidean")
print(fviz_cluster(clusters1, data = tabla.scaled, 
                   ellipse.type = "norm",ggtheme = theme_minimal())+
        labs(title = "Agrupamiento en 2 Clústers con distancia Euclideana"))

# Comprobando el mejor número de clusters por la distancia "manhattan" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
posibles.clusters2 <- NbClust(tabla.scaled, distance = "manhattan", min.nc=2, max.nc=10, method="kmeans",index="alllong")
print(fviz_nbclust(posibles.clusters2)+
        labs(title = "K óptimo con distancia Manhattan"))

# El mejor candidato es de 2 clusters, por tanto, se imprimen sus gráficos.
clusters2 <- pam(tabla.scaled, k = 2, metric = "manhattan")
print(fviz_cluster(clusters2, data = tabla.scaled, 
                   ellipse.type = "norm",ggtheme = theme_minimal())+
        labs(title = "Agrupamiento en 2 Clústers con distancia Manhattan"))

# Comparando ambos clusters, la mejor opción es el obtenido por la distancia de
# "manhattan", dado que este tiene menos datos mezclados entre sus grupos.

# Comprobando el mejor número de clusters por la distancia "gower" y el método
# "kmeans", entregando todos los índices posibles y ver la mejor opción.
cluster.daisy <- daisy(tabla.scaled, metric = "gower")
cluster.daisy.matrix <- as.matrix(cluster.daisy)
clusters3 <- kmeans(cluster.daisy.matrix, 2)
pam_fit <- pam(cluster.daisy, diss = TRUE, k = 2)

print(fviz_cluster(clusters3, data = tabla.scaled,
                   ellipse.type = "norm",ggtheme = theme_minimal())+
        labs(title = "Agrupamiento en 2 Clústers con distancia Gower"))

# Comparando los 3 clusters, la mejor opción es el obtenido por la distancia de
# "gower", dado que este tiene menos datos mezclados entre sus grupos. El 
# segundo mejor cluster el el obtenido por la distancia manhattan.

#-------------------------------------------------------------------------------

dist.eucl = dist(tabla.scaled, method = "euclidean")
print(dist.eucl)
fviz_dist(dist.eucl)

dist.manh = dist(tabla.scaled, method = "manhattan")
print(dist.manh)
print(fviz_dist(dist.manh))

gower.dist <- daisy(tabla.scaled, metric = "gower", stand = FALSE)
print(gower.dist)
print(fviz_dist(gower.dist))

gower.matrix <- as.matrix(gower.dist)

sil_width <- c(NA)
for(i in 2:10){  
 pam_fit <- pam(gower.dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}
print(plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width"))
lines(1:10, sil_width)

#Se trabjará con 2 cluster, ya que se sabe que son 2 valores en clases (comestible y venenoso)
#Entonces se quiere ver, cuantos valores pertenecen a cada grupo (cluster)
k<-2


#-------------------------------------------------------------------------------

tabla$clusters <- clusters3$cluster
print(summary(tabla[tabla["clusters"] == 1, ]))
print(summary(tabla[tabla["clusters"] == 2, ]))
print(summary(tabla["clusters"]))

# Se obtienen los medioides reales respecto a cada variable por cluster.
# Se crea una tabla que muestra los medioides de cada variables según los grupos
# de cluster por "gower".
pam_results <- tabla %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))

med.gower <- tabla[pam_fit$medoids,]
med.gower <- select(med.gower, -clusters)
medoids.1 <- med.gower[1,1]
medoids.2 <- med.gower[2,1]
cont <- 2
while (length(medoids.1)<10){
  medoids.1 <- c(medoids.1, med.gower[1,cont])
  medoids.2 <- c(medoids.2, med.gower[2,cont])
  cont <- cont + 1
}

medoids.1[1] <- "no-recurrence-events"
medoids.2[1] <- "no-recurrence-events"

medoids.1[2] <- "40-39"
medoids.2[2] <- "40-39"

medoids.1[3] <- "ge40"
medoids.2[3] <- "ge40"

medoids.1[4] <- "25-29"
medoids.2[4] <- "20-24"

medoids.1[5] <- "0-2"
medoids.2[5] <- "0-2"

medoids.1[6] <- "no"
medoids.2[6] <- "no"

medoids.1[8] <- "right"
medoids.2[8] <- "left"

medoids.1[9] <- "left_up"
medoids.2[9] <- "left_low"

medoids.1[10] <- "no"
medoids.2[10] <- "no"
# Se genera la tabla de medioides.
g.medoids <- data.frame(variables = columns,medoids.1, medoids.2)
print(g.medoids)
