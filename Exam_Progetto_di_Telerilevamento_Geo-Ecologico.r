################ MONITORAGGIO GEO-ECOLOGICO DEL BASSO SALENTO ##################
############## PROGETTO DI TELERILEVAMENTO GEO-ECOLOGICO - 93457 ###############



# Installazione dei packages necessari
install.packages("raster")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("viridis")

# Caricamento dei packages installati
library(raster)
library(ggplot2)
library(patchwork)
library(viridis)

# Impostare la working directory
setwd("C:/Data_telerilevamento/lab/Exam_Sentinel")





############# 1. IMPORTAZIONE E VISUALIZZAZIONE IMMAGINI Sentinel2 #############



### IMMAGINI 2018

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2018 <- list.files(pattern = "T34TBK_20180714T094029_B")
rlist_2018

# Applico la funzione raster() all'intera lista
import_2018 <- lapply(rlist_2018, raster)
import_2018

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2018 <- stack(import_2018)

# Visualizzo le informazioni
img_2018

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2018
plot(img_2018)

# Ritaglio l'area di interesse
ext <- c(246000,280000,4430000,4464100)
Sal_2018 <- crop(img_2018, ext)

# Plot di Sal_2018: colori reali e NIR + esportazione in .pdf
pdf("Salento_2018.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2018,4,3,2, stretch = "lin", main = "2018_TC")
plotRGB(Sal_2018,5,4,3, stretch = "lin", main = "2018_NIR")
dev.off()



### IMMAGINI 2019

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2019 <- list.files(pattern = "T34TBK_20190709T094039_B")
rlist_2019

# Applico la funzione raster() all'intera lista
import_2019 <- lapply(rlist_2019, raster)
import_2019

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2019 <- stack(import_2019)

# Visualizzo le informazioni
img_2019

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2019
plot(img_2019)

# Ritaglio l'area di interesse
Sal_2019 <- crop(img_2019, ext)

# Plot di Sal_2019: colori reali e NIR + esportazione in .pdf
pdf("Salento_2019.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2019,4,3,2, stretch = "lin", main = "2019_TC")
plotRGB(Sal_2019,5,4,3, stretch = "lin", main = "2019_NIR")
dev.off()



### IMMAGINI 2020

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2020 <- list.files(pattern = "T34TBK_20200723T094039_B")
rlist_2020

# Applico la funzione raster() all'intera lista
import_2020 <- lapply(rlist_2020, raster)
import_2020

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2020 <- stack(import_2020)

# Visualizzo le informazioni
img_2020

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2020
plot(img_2020)

# Ritaglio l'area di interesse
Sal_2020 <- crop(img_2020, ext)

# Plot di Sal_2020: colori reali e NIR + esportazione in .pdf
pdf("Salento_2020.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2020,4,3,2, stretch = "lin", main = "2020_TC")
plotRGB(Sal_2020,5,4,3, stretch = "lin", main = "2020_NIR")
dev.off()



### IMMAGINI 2021

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2021 <- list.files(pattern = "T34TBK_20210728T094029_B")
rlist_2021

# Applico la funzione raster() all'intera lista
import_2021 <- lapply(rlist_2021, raster)
import_2021

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2021 <- stack(import_2021)

# Visualizzo le informazioni
img_2021

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2021
plot(img_2021)

# Ritaglio l'area di interesse
Sal_2021 <- crop(img_2021, ext)

# Plot di Sal_2021: colori reali e NIR + esportazione in .pdf
pdf("Salento_2021.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2021,4,3,2, stretch = "lin", main = "2021_TC")
plotRGB(Sal_2021,5,4,3, stretch = "lin", main = "2021_NIR")
dev.off()



### IMMAGINI 2022

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2022 <- list.files(pattern = "T34TBK_20220703T094039_B")
rlist_2022

# Applico la funzione raster() all'intera lista
import_2022 <- lapply(rlist_2022, raster)

import_2022

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2022 <- stack(import_2022)

# Visualizzo le informazioni
img_2022

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2022
plot(img_2022)

# Ritaglio l'area di interesse
Sal_2022 <- crop(img_2022, ext)

# Plot di Sal_2022: colori reali e NIR + esportazione in .pdf
pdf("Salento_2022.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2022,4,3,2, stretch = "lin", main = "2022_TC")
plotRGB(Sal_2022,5,4,3, stretch = "lin", main = "2022_NIR")
dev.off()



### IMMAGINI 2023

# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2023 <- list.files(pattern = "T34TBK_20230728T093549_B")
rlist_2023

# Applico la funzione raster() all'intera lista
import_2023 <- lapply(rlist_2023, raster)
import_2023

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2023 <- stack(import_2023)

# Visualizzo le informazioni
img_2023

# Banda 2  = blu
# Banda 3  = verde
# Banda 4  = rosso
# Banda 8a = NIR

# Plot img_2023
plot(img_2023)

# Ritaglio l'area di interesse
Sal_2023 <- crop(img_2023, ext)

# Plot di Sal_2023: colori reali e NIR + esportazione in .pdf
pdf("Salento_2023.pdf")
par(mfrow = c(1,2))
plotRGB(Sal_2023,4,3,2, stretch = "lin", main = "2023_TC")
plotRGB(Sal_2023,5,4,3, stretch = "lin", main = "2023_NIR")
dev.off()



### PLOT DELLE 6 IMMAGINI + ESPORTAZIONE IN .pdf

# Colori reali
pdf("Salento_TC.pdf")
par(mfrow = c(2,3))
plotRGB(Sal_2018,4,3,2, stretch = "lin", main = "2018_TC")
plotRGB(Sal_2019,4,3,2, stretch = "lin", main = "2019_TC")
plotRGB(Sal_2020,4,3,2, stretch = "lin", main = "2020_TC")
plotRGB(Sal_2021,4,3,2, stretch = "lin", main = "2021_TC")
plotRGB(Sal_2022,4,3,2, stretch = "lin", main = "2022_TC")
plotRGB(Sal_2023,4,3,2, stretch = "lin", main = "2023_TC")
dev.off()

# Colori NIR
pdf("Salento_NIR.pdf")
par(mfrow = c(2,3))
plotRGB(Sal_2018,5,4,3, stretch = "lin", main = "2018_NIR")
plotRGB(Sal_2019,5,4,3, stretch = "lin", main = "2019_NIR")
plotRGB(Sal_2020,5,4,3, stretch = "lin", main = "2020_NIR")
plotRGB(Sal_2021,5,4,3, stretch = "lin", main = "2021_NIR")
plotRGB(Sal_2022,5,4,3, stretch = "lin", main = "2022_NIR")
plotRGB(Sal_2023,5,4,3, stretch = "lin", main = "2023_NIR")
dev.off()





################### 2. CALCOLO E PLOT DEGLI INDICI SPETTRALI ###################



### NDVI (NORMALIZED DIFFERENCE VEGETATION INDEX)

# Calcolo del DVI (Difference Vegetation Index)
# DVI = NIR - rosso
DVI_2018 <- Sal_2018[[5]] - Sal_2018[[4]]
DVI_2019 <- Sal_2019[[5]] - Sal_2019[[4]]
DVI_2020 <- Sal_2020[[5]] - Sal_2020[[4]]
DVI_2021 <- Sal_2021[[5]] - Sal_2021[[4]]
DVI_2022 <- Sal_2022[[5]] - Sal_2022[[4]]
DVI_2023 <- Sal_2023[[5]] - Sal_2023[[4]]

# Creazione di una palette di colori
cl <- colorRampPalette(c("blue","darkgrey","yellow"))(100)

# Plot DVI
par(mfrow = c(2,3))
plot(DVI_2018, col = cl, main = "DVI_2018")
plot(DVI_2019, col = cl, main = "DVI_2019")
plot(DVI_2020, col = cl, main = "DVI_2020")
plot(DVI_2021, col = cl, main = "DVI_2021")
plot(DVI_2022, col = cl, main = "DVI_2022")
plot(DVI_2023, col = cl, main = "DVI_2023")
dev.off()

# Calcolo del NDVI (Normalized Difference Vegetation Index)
# NDVI = (NIR - rosso) / (NIR + rosso) = DVI / (NIR + rosso)
NDVI_2018 <- DVI_2018 / (Sal_2018[[5]] + Sal_2018[[4]])
NDVI_2019 <- DVI_2019 / (Sal_2019[[5]] + Sal_2019[[4]])
NDVI_2020 <- DVI_2020 / (Sal_2020[[5]] + Sal_2020[[4]])
NDVI_2021 <- DVI_2021 / (Sal_2021[[5]] + Sal_2021[[4]])
NDVI_2022 <- DVI_2022 / (Sal_2022[[5]] + Sal_2022[[4]])
NDVI_2023 <- DVI_2023 / (Sal_2023[[5]] + Sal_2023[[4]])

# Plot NDVI
par(mfrow = c(2,3))
plot(NDVI_2018, col = cl, main = "NDVI_2018")
plot(NDVI_2019, col = cl, main = "NDVI_2019")
plot(NDVI_2020, col = cl, main = "NDVI_2020")
plot(NDVI_2021, col = cl, main = "NDVI_2021")
plot(NDVI_2022, col = cl, main = "NDVI_2022")
plot(NDVI_2023, col = cl, main = "NDVI_2023")
dev.off()

# Calcolo la differenza fra NDVI_2018 e NDVI_2023
NDVI_def <- NDVI_2018 - NDVI_2023

# Plot di NDVI_def + esportazione in .pdf
pdf("NDVI_def.pdf")
plot(NDVI_def, col = cl, main = "NDVI_def")
dev.off()

# GRIGIO: aree con differenza di NDVI fra 2018 e 2023 nulla
# GIALLO: perdita in copertura vegetale dal 2018 al 2023
# BLU:    guadagno in copertura vegetale dal 2018 al 2023



### MSAVI (MODIFIED SOIL-ADJUSTED VEGETATION INDEX)

# Indice introdotto per ridurre l'effetto del suolo nelle misurazioni di vegeta-
# zione. Rispetto all'NDVI, il MSAVI tende a essere meno sensibile all'effetto 
# del suolo e può fornire stime più accurate della copertura vegetale.

# Calcolo del MSAVI (Modified Soil-Adjusted Vegetation Index)
# MSAVI = {2*NIR + 1 - sqrt[(2*NIR + 1)^2 - 8*(NIR - rosso)]} / 2
MSAVI_2018 <- ((2*Sal_2018[[5]] + 1) - sqrt((2*Sal_2018[[5]] + 1)^2 - 8 *(Sal_2018[[5]] - Sal_2018[[4]]))) / 2
MSAVI_2019 <- ((2*Sal_2019[[5]] + 1) - sqrt((2*Sal_2019[[5]] + 1)^2 - 8 *(Sal_2019[[5]] - Sal_2019[[4]]))) / 2
MSAVI_2020 <- ((2*Sal_2020[[5]] + 1) - sqrt((2*Sal_2020[[5]] + 1)^2 - 8 *(Sal_2020[[5]] - Sal_2020[[4]]))) / 2
MSAVI_2021 <- ((2*Sal_2021[[5]] + 1) - sqrt((2*Sal_2021[[5]] + 1)^2 - 8 *(Sal_2021[[5]] - Sal_2021[[4]]))) / 2
MSAVI_2022 <- ((2*Sal_2022[[5]] + 1) - sqrt((2*Sal_2022[[5]] + 1)^2 - 8 *(Sal_2022[[5]] - Sal_2022[[4]]))) / 2
MSAVI_2023 <- ((2*Sal_2023[[5]] + 1) - sqrt((2*Sal_2023[[5]] + 1)^2 - 8 *(Sal_2023[[5]] - Sal_2023[[4]]))) / 2

# Plot MSAVI
par(mfrow = c(2,3))
plot(MSAVI_2018, col = cl, main = "MSAVI_2018")
plot(MSAVI_2019, col = cl, main = "MSAVI_2019")
plot(MSAVI_2020, col = cl, main = "MSAVI_2020")
plot(MSAVI_2021, col = cl, main = "MSAVI_2021")
plot(MSAVI_2022, col = cl, main = "MSAVI_2022")
plot(MSAVI_2023, col = cl, main = "MSAVI_2023")
dev.off()

# Calcolo la differenza fra MSAVI_2018 e MSAVI_2023
MSAVI_def <- MSAVI_2018 - MSAVI_2023

# Plot di NDVI_def + esportazione in .pdf
pdf("MSAVI_def.pdf")
plot(MSAVI_def, col = cl, main = "MSAVI_def")
dev.off()

# GRIGIO: aree con differenza di NDVI fra 2018 e 2023 nulla
# GIALLO: perdita in copertura vegetale dal 2018 al 2023
# BLU:    guadagno in copertura vegetale dal 2018 al 2023



### VERIFICO SE NDVI_def e MSAVI_def DIFFERISCONO SIGNIFICATIVAMENTE

# Plot di NDVI_def e MSAVI_def + esportazione in .pdf
pdf("NDVI_def + MSAVI_def.pdf")
par(mfrow = c(1,2))
plot(NDVI_def, col = cl, main = "NDVI_def")
plot(MSAVI_def, col = cl, main = "MSAVI_def")
dev.off()

# Eseguo un T-test su NDVI_def e MSAVI_def
t_result <- t.test(NDVI_def[], MSAVI_def[], paired = TRUE)
t_result





#################### 3. PCA (PRINCIPAL COMPONENT ANALYSIS) #####################



# ImpostO il seme del generatore di numeri casuali
set.seed(1)

# Unisco gli indici NDVI_def e MSAVI_def in un unico oggetto
box <- stack(NDVI_def, MSAVI_def)

# Plot
plot(box, main = "Differenze fra NDVI_def e MSAVI_def", xaxt = "n", yaxt = "n")

# Effettuo un campionamento casuale di 10000 pixel da box
sr <- sampleRandom(box, 10000)

# Effettuo la PCA (Principal Component Analysis)
PCA <- prcomp(sr)

# Visualizzazione delle informazioni relative alla PCA
summary(PCA)

# Plot della varianza spiegata da ciascuna delle componenti
plot(PCA)

# Proiezione dell'oggetto box nello spazio creato precedentemente usando le CP
PCI <- predict(box, PCA, index = 1:2)

# Plot della PC1
plot(PCI[[1]])

# Conversione di PC1 in un dataframe
PC_fin <- as.data.frame(PCI[[1]], xy = T)

# Plot 
plot_PCA <- ggplot() + 
            geom_raster(PC_fin, mapping = aes(x = x, y = y, fill = PC1)) + 
            scale_fill_viridis(option="magma") +
            labs(title = "PC1")

# Maggior variabilità = valore PC1 più basso
# Minor variabilità   = valore PC1 più alto

# Esportazione di final_plot in .pdf
pdf("PC1.pdf")
print(plot_PCA)
dev.off()





################################ 4. LAND COVER #################################



### CLASSIFICAZIONE IMMAGINI 2018

# Estrazione dei valori dalle immagini del 2018
single_nr_2018 <- getValues(Sal_2018)
single_nr_2018

# Classificazione
k_cluster_2018 <- kmeans(single_nr_2018, centers=2)
k_cluster_2018

# Set dei valori
Sal_2018_class <- setValues(Sal_2018[[1]], k_cluster_2018$cluster)

# Creazione di una palette di colori
cl_freq <- colorRampPalette(c("blue","yellow"))(2)

# Plot
plot(Sal_2018_class, col = cl_freq)

# Calcolo delle frequenze
freq_2018 <- freq(Sal_2018_class)
freq_2018

# Calcolo pixel totali
tot <- ncell(Sal_2018_class)
tot

# Calcolo delle percentuali
perc_2018 <- round((freq_2018 * 100) / tot, digit=5)
perc_2018



### CLASSIFICAZIONE IMMAGINI 2019

# Estrazione dei valori dalle immagini del 2019
single_nr_2019 <- getValues(Sal_2019)
single_nr_2019

# Classificazione
k_cluster_2019 <- kmeans(single_nr_2019, centers=2)
k_cluster_2019

# Set dei valori
Sal_2019_class <- setValues(Sal_2019[[1]], k_cluster_2019$cluster)

# Plot
plot(Sal_2019_class, col = cl_freq)

# Calcolo delle frequenze
freq_2019 <- freq(Sal_2019_class)
freq_2019

# Calcolo delle percentuali
perc_2019 <- round((freq_2019 * 100) / tot, digit=5)
perc_2019



### CLASSIFICAZIONE IMMAGINI 2020

# Estrazione dei valori dalle immagini del 2020
single_nr_2020 <- getValues(Sal_2020)
single_nr_2020

# Classificazione
k_cluster_2020 <- kmeans(single_nr_2020, centers=2)
k_cluster_2020

# Set dei valori
Sal_2020_class <- setValues(Sal_2020[[1]], k_cluster_2020$cluster)

# Plot
plot(Sal_2020_class, col = cl_freq)

# Calcolo delle frequenze
freq_2020 <- freq(Sal_2020_class)
freq_2020

# Calcolo delle percentuali
perc_2020 <- round((freq_2020 * 100) / tot, digit=5)
perc_2020



### CLASSIFICAZIONE IMMAGINI 2021

# Estrazione dei valori dalle immagini del 2021
single_nr_2021 <- getValues(Sal_2021)
single_nr_2021

# Classificazione
k_cluster_2021 <- kmeans(single_nr_2021, centers=2)
k_cluster_2021

# Set dei valori
Sal_2021_class <- setValues(Sal_2021[[1]], k_cluster_2021$cluster)

# Plot
plot(Sal_2021_class, col = cl_freq)

# Calcolo delle frequenze
freq_2021 <- freq(Sal_2021_class)
freq_2021

# Calcolo delle percentuali
perc_2021 <- round((freq_2021 * 100) / tot, digit=5)
perc_2021



### CLASSIFICAZIONE IMMAGINI 2022

# Estrazione dei valori dalle immagini del 2022
single_nr_2022 <- getValues(Sal_2022)
single_nr_2022

# Classificazione
k_cluster_2022 <- kmeans(single_nr_2022, centers=2)
k_cluster_2022

# Set dei valori
Sal_2022_class <- setValues(Sal_2022[[1]], k_cluster_2022$cluster)

# Plot
plot(Sal_2022_class, col = cl_freq)

# Calcolo delle frequenze
freq_2022 <- freq(Sal_2022_class)
freq_2022

# Calcolo delle percentuali
perc_2022 <- round((freq_2022 * 100) / tot, digit=5)
perc_2022



### CLASSIFICAZIONE IMMAGINI 2023

# Estrazione dei valori dalle immagini del 2023
single_nr_2023 <- getValues(Sal_2023)
single_nr_2023

# Classificazione
k_cluster_2023 <- kmeans(single_nr_2023, centers=2)
k_cluster_2023

# Set dei valori
Sal_2023_class <- setValues(Sal_2023[[1]], k_cluster_2023$cluster)

# Plot
plot(Sal_2023_class, col = cl_freq)

# Calcolo delle frequenze
freq_2023 <- freq(Sal_2023_class)
freq_2023

# Calcolo delle percentuali
perc_2023 <- round((freq_2023 * 100) / tot, digit=5)
perc_2023



### VISUALIZZAZIONE DEI RISULTATI


perc_2018
#        value    count
# [1,] 0.00031 82.28862
# [2,] 0.00062 17.71138

perc_2019
#        value    count
# [1,] 0.00031 73.49073
# [2,] 0.00062 26.50927

perc_2020
#        value    count
# [1,] 0.00031 74.11995
# [2,] 0.00062 25.88005

perc_2021
#        value    count
# [1,] 0.00031 70.56193
# [2,] 0.00062 29.43807

perc_2022
#        value    count
# [1,] 0.00031 68.82305
# [2,] 0.00062 31.17695

perc_2023
#        value    count
# [1,] 0.00031 67.17314
# [2,] 0.00062 32.82686



### CREAZIONE DI UN DATAFRAME CON I RISULTATI DELLA CLASSIFICAZIONE

# Creazione dei vettori
copertura_vegetale <- c("buona","ridotta/assente")
P_2018 <- c(82.28862, 17.71138)
P_2019 <- c(73.49073, 26.50927)
P_2020 <- c(74.11995, 25.88005)
P_2021 <- c(70.56193, 29.43807)
P_2022 <- c(68.82305, 31.17695)
P_2023 <- c(67.17314, 32.82686)

# Creazione del dataframe
Land_cover_perc <- data.frame(copertura_vegetale, P_2018, P_2019, P_2020, P_2021, P_2022, P_2023)

# Visualizzazione del dataframe
Land_cover_perc



### PLOTE DELLE PERCENTUALI DI LAND COVER + ESPORTAZIONE IN .pdf

cl_barplot <- c("buona" = "blue", "ridotta/assente" = "yellow")

# 2018
pdf("Percentuali Land Cover 2018.pdf")
plot_2018 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale,y = P_2018, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) +
             labs(x = "Land Cover", y = "%", title = "Land Cover 2018") +
             theme(legend.position = "none") +
             ylim(0, 100) +
             geom_text(aes(label = sprintf("%.2f%%", P_2018), y = P_2018), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# 2019
pdf("Percentuali Land Cover 2019.pdf")
plot_2019 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale,y = P_2019, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) +
             labs(x = "Land Cover", y = "%", title = "Land Cover 2019") +
             theme(legend.position = "none") +
             ylim(0, 100) +
             geom_text(aes(label = sprintf("%.2f%%", P_2019), y = P_2019), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# 2020
pdf("Percentuali Land Cover 2020.pdf")
plot_2020 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale,y = P_2020, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) +
             labs(x = "Land Cover", y = "%", title = "Land Cover 2020") +
             theme(legend.position = "none") +
             ylim(0, 100) +
             geom_text(aes(label = sprintf("%.2f%%", P_2020), y = P_2020), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# 2021
pdf("Percentuali Land Cover 2021.pdf")
plot_2021 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale, y = P_2021, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) +
             labs(x = "Land Cover", y = "%", title = "Land Cover 2021") +
             theme(legend.position = "none") +
             ylim(0, 100) + 
             geom_text(aes(label = sprintf("%.2f%%", P_2021), y = P_2021), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# 2022
pdf("Percentuali Land Cover 2022.pdf")
plot_2022 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale, y = P_2022, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) +
             labs(x = "Land Cover", y = "%", title = "Land Cover 2022") +
             theme(legend.position = "none") +
             ylim(0, 100) + 
             geom_text(aes(label = sprintf("%.2f%%", P_2022), y = P_2022), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# 2023
pdf("Percentuali Land Cover 2023.pdf")
plot_2023 <- ggplot(Land_cover_perc,aes(x = copertura_vegetale, y = P_2023, fill = copertura_vegetale)) +
             geom_bar(stat = "identity") +
             scale_fill_manual(values = cl_barplot) + 
             labs(x = "Land Cover", y = "%", title = "Land Cover 2023") +
             theme(legend.position = "none") +
             ylim(0, 100) + 
             geom_text(aes(label = sprintf("%.2f%%", P_2023), y = P_2023), 
             position = position_stack(vjust = 0.5), size = 4)
dev.off()

# Unione dei quattro plot in un unico grafico
final_plot <- plot_2018 + plot_2019 + plot_2020 + plot_2021 + plot_2022 + plot_2023
final_plot

# Esportazione di final_plot in .pdf
pdf("Final Plot.pdf")
print(final_plot)
dev.off()

