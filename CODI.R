# CODI NO INCLÓS QUE PODRIA SER ÚTIL. 

getwd()
wd<-"/Users/annia/Library/CloudStorage/GoogleDrive-acniell@gmail.com/My Drive/UOC/TERCER QUADRIMESTRE/ANÀLISI DADES ÒMIQUES/PEC2/PEC2"
dataDir <- file.path(wd, "dades")
resultsDir <- file.path(wd, "results")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("affy", "Biobase"))
BiocManager::install("GEOquery")
if (!requireNamespace("oligo", quietly = TRUE)) {
  BiocManager::install("oligo")
}
library(GEOquery)
library(oligo)
library(affy)
library(Biobase)

GSE38531<-getGEO("GSE38531", GSEMatrix = TRUE)
GSE38531

allTargets <- read.table("allTargets.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

filter_microarray <- function(allTargets, seed =123) {
  set.seed(seed)
  filtered <- subset(allTargets, time != "hour 2")
  
  # Dividir el dataset por grupos únicos de 'infection' + 'agent'
  filtered$group <- interaction(filtered$infection, filtered$agent)
  
  # Seleccionar 4 muestras al azar de cada grupo
  selected <- do.call(rbind, lapply(split(filtered, filtered$group), function(group_data) {
    if (nrow(group_data) > 4) {
      group_data[sample(1:nrow(group_data), 4), ]
    } else {
      group_data
    }
  }))
  
  # Obtener los índices originales como nombres de las filas seleccionadas
  original_indices <- match(selected$sample, allTargets$sample)
  
  # Modificar los rownames usando 'sample' y los índices originales
  rownames(selected) <- paste0(selected$sample, ".", original_indices)
  
  # Eliminar la columna 'group' y devolver el resultado
  selected$group <- NULL
  return(selected)
}

result <- filter_microarray(allTargets, seed=41577436)
print(result)
write.table(result, "mostra_seleccionada_41577436.txt", sep = "\t", row.names = FALSE)

GSE38531<-getGEO("GSE38531", GSEMatrix = TRUE)
GSE38531 #L'expressionSet

ES<- GSE38531[[1]] 
matriu_expressio<-exprs(ES)
head(matriu_expressio)
fenodata<-pData(ES)
head(fenodata)

mostres <- read.table("mostra_seleccionada_41577436.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
noms_mostra<-mostres$sample #llista dels que volem. 

sampleNames(ES)->noms_ES
coincidencies<-intersect(noms_mostra, noms_ES)
coincidencies #veiem que efectivament hi ha 24 elements a la llista. 
ES2<-ES[,coincidencies]
ES2


