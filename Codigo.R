#Trabajo final
#instalar devtools y paquetes necesarios
library(devtools)
library(tidyverse)
library(dplyr)

# Primer Utiliza la función merge para combinar las matrices BRCA_normal y BRCA_PT
# basándose en las columnas en comun, al final pido que  incluir todas las filas,
# tanto las que tienen coincidencias como las que no.
datos_concardenados <- merge(DynamicCancerDriverKM::BRCA_normal,
                             DynamicCancerDriverKM::BRCA_PT,all = TRUE)

#Resumo para entender un poco mejor los datos pero aun no es muy concluyente
summary(datos_concardenados)
#intento sacar la media pero tengo variables categoricas
median(datos_concardenados)

# Es por esto que elimino las columas que no me interesan y dejo solo genes y sample
datos_concardenados <- datos_concardenados[, -c(1:3, 5:7)]

## Tenia que sacar el numero mas alto entre los datos, exclui a la primera columna (sample)
# y multiplico por 0.00035 que es el 0.035% , na para excluir estos valores, con esto
# puedo revisar que genes estan expresados
#calculo la media para tener un punto de partida y luego sacar un valor, lo intente solo
#con la media aritmetica pero me enrrede, por eso tuve que calcular un valor y multiplicar
media_datos <- mean(unlist(datos_concardenados[, -1]))
# escogi el 0.0003 por estar cerca a la media
numero_mayor <- max(datos_concardenados[,2:23688], na.rm=FALSE)*0.0003


### Creo una copia donde pondre 1 a los datos que pasen el umbral establecido y cero los que no

binarizacion_datos <- datos_concardenados

# utilice apply para que tome cada columna x y utilice ifelse para asignar 1 si el valor
# es mayor queel numero establecido como umbral

binarizacion_datos[, 2:23688] <- apply(datos_concardenados[, 2:23688], 2,
                                       function(x) ifelse(x > numero_mayor, 1, 0))

## Limpie los datos  el eliminar losgenes que no estan expresados por lo menos en el 20%

no_expresados <- colMeans(binarizacion_datos[, 2:ncol(binarizacion_datos)] == 1)

Variables_expresadas <- binarizacion_datos[, no_expresados >= 0.2]

# Calcule las medias de cada columna de la matriz booleana resultante, donde la media
# representa el porcentaje de unos en cada columna.
medias_variables <- colMeans(Variables_expresadas[, 2:ncol(Variables_expresadas)]==1)
medias_variables

# Cargue PPI que contiene todos los genes
PPI <- DynamicCancerDriverKM::PPI
# renombre las variables para que concuerden respecto a HGNC.symbol para poder concardenarlas
titulos_ppi <- AMCBGeneUtils::changeGeneId(PPI$`Input-node Gene Symbol`,from="HGNC.symbol")

# nuevamente para acceder por columnas con colnames
colnames(Variables_expresadas)[c(2:ncol(Variables_expresadas))] <- titulos_ppi$HGNC.symbol

# Con la funcion group_by agrupo la cantidad de veces que aparece el nombre de un gen, comi
# un filtro de excel y les asigno el mismo nombre para luego agrupar en un resultante

Grupo_in<- PPI %>% group_by(`Input-node Gene Symbol`) %>%
  summarise(NN = n()) %>%
  rename(`PPI_genes`=`Input-node Gene Symbol`)

Grupo_out <- PPI %>% group_by(`Output-node Gene Symbol`) %>%
  summarise(NN = n()) %>%
  rename(`PPI_genes`=`Output-node Gene Symbol`)

### Agrupo en un resultante mediante la columna PPI_genes
Resultantes <- bind_rows(Grupo_in, Grupo_out) %>%
  group_by(`PPI_genes`) %>%
  summarise(NN = sum(NN)) %>%
  arrange(desc(NN)) %>% #organizo en orden desc
  top_n(100, NN) #Saco el top 100 de las variables con mayor presencia

#CAMBIO LOS TITULOS EN FUNCION DE HGNC
titulos_resultantes <- AMCBGeneUtils::changeGeneId(Resultantes$PPI_genes,from="HGNC.symbol")
rownames(Resultantes) <- titulos_resultantes$HGNC.symbol

#Accedo a los titulos de cada variable genica

Cambio_titulos <- colnames(Variables_expresadas)
# Convertir los nombres de las variables en una matriz para unir con las resultantes
Cambio_titulos_matriz <- as.matrix(Cambio_titulos)

interseccion <- intersect(Resultantes$PPI_genes, Cambio_titulos_matriz)
interseccion <- as.matrix(interseccion)

#####
Matriz_interseccion <- Variables_expresadas[,c("sample_type",interseccion)]


#########################################################################
###################### MODELO KNN ########################################
#########################################################################
install.packages("caret")
library(caret)

# duplico primero mis datos
interseccion_datos <- Matriz_interseccion
# Ahora convierto sample_type en un factor
interseccion_datos$sample_type <- as.factor(interseccion_datos$sample_type)

set.seed(1)# Establezco una semilla para garantizar la reproducibilidad de los resultados
# cuando se generan números aleatorios. Cada que ejecute el mismo código, obtengo los
# mismos resultados.


#mediante muestreo aleatorio tomo un 70% para entrenamiento
sample.index <- sample(1:nrow(interseccion_datos),nrow(interseccion_datos)*0.7
                       ,replace = F)

#creo una variable para guardar los titulos de las columnas menos el de sample_type
predictor <- colnames(interseccion_datos)[-1]

#creo el conjunto de datos con el 70% que tome del muestreo aleatorio
train.data <- interseccion_datos[sample.index,c(predictor,"sample_type"),drop=F]
#creo un conjuntos de datos con el 30% restante
test.data <- interseccion_datos[-sample.index,c(predictor,"sample_type"),drop=F]

# Establezco control para el entrenamiento del modelo utilizando validación cruzada
# con una proporción del 70% para entrenamiento.
ctrl <- trainControl(method="cv", p = 0.7)

# Entreno el modelo utilizando los datos de entrenamiento. Se indica que se realizará
# una validación cruzada y se aplicará procesamiento previo como centrado y escalado.
# Se busca la longitud óptima del vecino (k) en el rango de 1 a 25.

knn_train <- train( sample_type ~ .
                   , data = train.data
                   , method = "knn", trControl = ctrl
                   , preProcess = c("center","scale")
                   , tuneLength = 25)

knn_train
plot(knn_train)

# Realiza predicciones en los datos de prueba utilizando el modelo entrenado.
knn_predictor <- predict(knn_train, newdata = test.data)
knn_predictor

confusionMatrix(knn_predictor, test.data$sample_type)










