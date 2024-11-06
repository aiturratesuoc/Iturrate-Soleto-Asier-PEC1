library(SummarizedExperiment)
library(ggfortify)
library(pheatmap)

## Importar dataset y adapatación

# Cargamos el dataset
human_cachexia <- read.csv("human_cachexia.csv", header = TRUE)

# Visualizamos datase
head(human_cachexia)

# Hacemos un resumen de la estructura de nuestro dataset
str(human_cachexia)

# Cambiamos los ID de los pacientes por paciente+numeros
human_cachexia$Patient.ID[1:47] <- paste0("paciente", 1:47)

# Cambiamos los ID de los controles por control+numero
human_cachexia$Patient.ID[(nrow(human_cachexia)-29):nrow(human_cachexia)] <- paste0("control", 1:30)



## Generamos matriz de conteos

# Generamos la matriz transpuesta de nuestro dataset con los datos de los metabolitos
myCounts <- data.frame(t(human_cachexia[3:65]))

# Añadimos el nombre de los individuos a cada columna
colnames(myCounts) <- human_cachexia[, 1]

# Visualizamos la matriz de conteos
head(myCounts)


## Generamos matriz de covariables

# Extraemos las dos primeras columnas de nuestro dataset original
myGroups <- data.frame(sampleName = human_cachexia[, 1], 
                     group = human_cachexia[, 2], 
                     row.names = 1)

# Visualizamos matriz de covariables
head(myGroups)


## Generamos vector con identificadores

# Guardamos los nombres de los metabolitos
myMetabolites <- c(colnames(human_cachexia[3:65]))


## Generamos información del SummaryzedExperiment
myInfo <- list(myName = "Asier Iturrate", 
               myInstitution = "UOC", 
               myContact = "aiturrates@uoc.edu", 
               myTitle = "SummarizedExperiment generado para la PEC1",
               Description = "Objeto original creado a partir del dataset human_cachexia" )



## Generamos SummaryzedExperiment

# Creamos SummarizedExperiment
mySE <- SummarizedExperiment(assays = list(metabolomica = myCounts),
                             colData = myGroups, 
                             rowData = myMetabolites)

# Añadimos metadatos
metadata(mySE) <- myInfo

# Visualizamos información de SummarizedExperiment
show(mySE)

# Guardamos el SummaryzedExperiment con formato .Rda
save(mySE, file = "mySE.Rda")


## Extraemos información del objeto SummaryzedExperiment

# Extraemos matriz de conteos
countData <- assays(mySE)$metabolomica

# Extraemos matriz de grupos
groupData <- as.data.frame(colData(mySE))



## Análisis exploratorio: PCA

# Calculamos los valores de cada componente principal
pca <- prcomp(t(countData), scale. = TRUE)

# Ploteamos PCA
autoplot(pca, 
         data = groupData,
         colour = "group",
         label = TRUE,
         label.size = 4,
         main = "PCA")



## Análisis exploratorio: clusterización jerárquica

# Normalizamos los datos por filas
scaled_countData <- scale(t(countData))  

# Realizamos el la clusterización jerárquica
hclust_result <- hclust(dist(scaled_countData))

# Ploteamos el dendograma
plot(hclust_result, 
     cex = 0.6,
     main = "Clusterización jerárquica",
     hang = -1)



## Análisis exploratorio: heatmap

pheatmap(countData, 
         scale = "row",
         annotation_col = groupData,
         cluster_rows = FALSE,
         cex = 0.8,
         main = "Heatmap")


## Actualizamos SummaryzedExperiment eliminando muestras

# Eliminamos las entradas de los pacientes correspondientes
mySE_new <- mySE[, !(colnames(mySE) %in% c("paciente27", "paciente4", "paciente5"))]


# Actualizamos las información del SummaryzedExperiment
myInfo <- list(myName = "Asier Iturrate", 
               myInstitution = "UOC", 
               myContact = "aiturrates@uoc.edu", 
               myTitle = "SummarizedExperiment filtrado generado para la PEC1",
               Description = "Se trata de un nuevo set de datos eliminando algunas entradas del dataset original")

# Actualizamos metadata
metadata(mySE_new) <- myInfo

# Guardamos el nuevo SummaryzedExperiment con formato .Rda
save(mySE_new, file = "mySE_new.Rda")


# Extraemos nueva matriz de conteos
countData_new <- assays(mySE_new)$metabolomica

# Extraemos nueva matriz de grupos
groupData_new <- as.data.frame(colData(mySE_new))


## Análisis exploratorio: PCA

# Calculamos los valores de cada componente principal
pca2 <- prcomp(t(countData_new), scale. = TRUE)

# Ploteamos PCA
autoplot(pca2, 
         data = groupData_new,
         colour = "group",
         label = TRUE,
         label.size = 4,
         main = "PCA nuevo dataset")



## Análisis exploratorio: clusterización jerárquica

# Normalizamos los nuevos datos por filas
scaled_countData_new <- scale(t(countData_new))  

# Realizamos  la clusterización jerárquica
hclust_result_new <- hclust(dist(scaled_countData_new))

# Ploteamos el dendograma
plot(hclust_result_new, 
     cex = 0.6,
     main = "Clusterización jerárquica nuevo dataset",
     hang = -1)



## Análisis exploratorio: heatmap
pheatmap(countData_new, 
         scale = "row",
         annotation_col = groupData_new,
         cluster_rows = FALSE,
         cex = 0.8,
         main = "Heatmap nuevo dataset")



