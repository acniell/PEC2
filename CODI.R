
knitr::opts_chunk$set(echo = TRUE)

install.packages("BiocManager")
BiocManager::install("oligo")
BiocManager::install("pd.mogene.2.1.st")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("pvca")
BiocManager::install("limma")
BiocManager::install("genefilter")
BiocManager::install("mogene21sttranscriptcluster.db")
BiocManager::install("annotate")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")
BiocManager::install(c("affy", "Biobase"))
BiocManager::install("GEOquery")
library(GEOquery)
library(oligo)
library(affy)
library(Biobase)
install.packages("knitr")
install.packages("colorspace")
install.packages("gplots")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("htmlTable")
install.packages("prettydoc")
install.packages("devtools")

getwd()
wd<-"/Users/annia/Library/CloudStorage/GoogleDrive-acniell@gmail.com/My Drive/UOC/TERCER QUADRIMESTRE/ANÀLISI DADES ÒMIQUES/PEC2/PEC2"
dataDir <- file.path(wd, "dades")
resultsDir <- file.path(wd, "results")

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

resultat_mostra <- filter_microarray(allTargets, seed=41577436)
print(resultat_mostra)

### Aprofitarem i ajustarem ara els noms, perquè no s'havia fet prèviament i ha generat problemes amb limma.

resultat_mostra$infection <- gsub("S\\. aureus USA300", "aureus", resultat_mostra$infection) 
resultat_mostra$time <- gsub("hour ", "", resultat_mostra$time)
resultat_mostra$agent <- gsub("linezolid", "line", gsub("vancomycin", "vanco", resultat_mostra$agent))  
resultat_mostra$nom <- paste0(
  sub("GSM944", "", resultat_mostra$sample), "_", 
  resultat_mostra$infection, "_", 
  resultat_mostra$time, "_", 
  resultat_mostra$agent
)

print(head(resultat_mostra))
write.table(resultat_mostra, file = file.path(dataDir, "mostra_seleccionada_41577436.txt"), 
            sep = "\t", row.names = FALSE)

GSE38531<-getGEO("GSE38531", GSEMatrix = TRUE)
GSE38531 #L'expressionSet

ES<- GSE38531[[1]] 
matriu_expressio<-exprs(ES)
head(matriu_expressio)
fenodata<-pData(ES)
head(fenodata)

mostres <- read.table(file.path(dataDir, "mostra_seleccionada_41577436.txt"), 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
noms_mostra<-mostres$sample #llista dels que volem. 

sampleNames(ES)->noms_ES
coincidencies<-intersect(noms_mostra, noms_ES)
coincidencies #veiem que efectivament hi ha 24 elements a la llista. 
ES2<-ES[,coincidencies]
ES2

#ens el guardem per si el necessitem més tard naive
write.csv(exprs(ES2), file = file.path(dataDir, "ES2.csv"))


library(oligo) #pels elements CET
library(Biobase) 
arxius_CEL <- list.celfiles("/Users/annia/Library/CloudStorage/GoogleDrive-acniell@gmail.com/My Drive/UOC/TERCER QUADRIMESTRE/ANÀLISI DADES ÒMIQUES/PEC2/PEC2/ARXIUS CEL", full.names = TRUE)
phenoData<-AnnotatedDataFrame(data = mostres)
rawData<-read.celfiles(arxius_CEL, phenoData=phenoData)
rawData

library(arrayQualityMetrics)
arrayQualityMetrics(rawData)

ESet<-rma(rawData)
ESet

sampleNames(ESet)<-mostres$sample
sampleNames(ESet)
pData(ESet)
head(ESet)

matriu2<-exprs(ESet)
dim(matriu2)
str(matriu2)
head(matriu2)
summary(matriu2)
sum(is.na(matriu2))

boxplot(matriu2, main="Boxplot", col = rainbow(ncol(matriu2)), cex.axis=0.4, las=2)
hist(as.numeric(as.matrix(matriu2)), 
     main="Representació dels valors d'expressió",
     xlab="Expressió gènica", 
     ylab="Freqüència d'expressió dels gens",
     col="cadetblue3", breaks=35)
distancia <- dist(t(matriu2)) 
clustering <- hclust(distancia)
plot(clustering, 
     main="Visualització de clústers",
     xlab="Mostres", ylab="Distància", 
     cex=0.5) 
pca <- prcomp(t(matriu2), scale.=TRUE) 
summary(pca)

plot(pca$x[,1], pca$x[,2], 
     xlab="PC1", 
     ylab="PC2",
     main="Anàlisi de components principals",
     pch=19, col="orange")

text(pca$x[,1], pca$x[,2], labels=colnames(matriu2), pos=1, cex=0.4)

desv_std<-apply(matriu2, 1, sd)
sorted_desv_std<-sort(desv_std)
head(sorted_desv_std)
plot(1:length(sorted_desv_std), sorted_desv_std, 
     main="Variabilitat dels gens",
     xlab="Índex de variabilitat", 
     ylab="SD",
     type="l", col="purple", lwd=1)

# Afegir línies verticals als percentils 90% i 95%
abline(v=length(sorted_desv_std)*c(0.9, 0.95), col=c("darkgreen", "cadetblue3"), lwd=1, lty=4)

# Llegenda
legend("topleft", legend=c("90%", "95%"), 
       col=c("darkgreen", "cadetblue3"), lwd=1, lty=4)

limit90<- quantile(desv_std, 0.9)
pool10 <- matriu2[desv_std >= limit90, ]
dim(pool10) 

mostres
str(mostres)
mostres$infection <- factor(mostres$infection, levels=c("uninfected", "aureus"))
mostres$time <- factor(as.character(mostres$time), levels=c("0", "24"))
mostres$agent <- factor(mostres$agent, levels=c("untreated", "line", "vanco"))
str(mostres)
mostres

disseny <- model.matrix(~ 0 + infection * time * agent, data = mostres)
colnames(disseny) <- gsub(":", "_", colnames(disseny))
disseny
colnames(disseny)

library(limma)
contrasts <- makeContrasts(
  Infectats_vs_Noinfectats_Untreated = infectionaureus - infectionuninfected,
  Infectats_vs_Noinfectats_Linezolid = infectionaureus_agentline - (infectionuninfected + agentline),
  Infectats_vs_Noinfectats_Vancomycin = infectionaureus_agentvanco - (infectionuninfected + agentvanco),
  levels = disseny
)

print(contrasts)
