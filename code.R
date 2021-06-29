library(KnowSeq)
library(CORElearn)
library(dplyr)
library(caret)
library(car)
library(R.utils)
require(ROCR)

# DATA READ AND PREPARATION

#MICROARRAY

# RAW microarray data is available in NCBI/GEO under ID mentioned in the main text. Once each serie was preprocess the integration between them
# was carried out applying normalizeBetweenArrays() limma function.

mic_exp <- read.csv('mic_exp.csv', header = TRUE)
mic_exp[,1] <- NULL
samples_mic <- read.csv('mic_samples.csv')
samples_mic[,1] <- NULL
expressionMatrix_tri_mic_1 <- as.matrix(mic_exp)
labels_tri_mic_1 <- samples_mic[,2]
expressionMatrix_ACC_mic_1 <- expressionMatrix_tri_mic_1[which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'ACC'),]
labels_ACC_mic_1 <- labels_tri_mic_1[which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'ACC')]
expressionMatrix_SCC_mic_1 <- expressionMatrix_tri_mic_1[which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'SCC'),]
labels_SCC_mic_1 <- labels_tri_mic_1[which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'SCC')]


# RNA-SEQ

#ACC

# RAW DATA PREPROCESS 
#countsInformation_ACC <- countsToMatrix('D:/documentos/master/TFM/Datos_counts/LUAD_pacientes.csv', extension = "")
#countsMatrix_ACC <- countsInformation_ACC$countsMatrix
#labels_ACC_rna <- countsInformation_ACC$labels 
#for(i in 1:length(labels_ACC_rna)){
  if (labels_ACC_rna[i] == 'Adenocarcinoma'){
    labels_ACC_rna[i] = 'ACC'
  }
  else {
    labels_ACC_rna[i] = 'Control'
  }
}
#myAnnotation_ACC <- getGenesAnnotation(rownames(countsMatrix_ACC))
#expressionMatrix_ACC_rna_1 <- calculateGeneExpressionValues(countsMatrix_ACC,myAnnotation_ACC,genesNames = TRUE)
#outliers_ACC_rna <- RNAseqQA(expressionMatrix_ACC_rna_1, toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE, outdir = '') 
#which(colnames(expressionMatrix_ACC_rna_1) %in% outliers_ACC_rna$outliers)
#expressionMatrix_ACC_rna <- expressionMatrix_ACC_rna[,-which(colnames(expressionMatrix_ACC_rna_1) %in% outliers_ACC_rna$outliers)]
#labels_ACC_rna <- labels_ACC_rna[-which(colnames(expressionMatrix_ACC_rna_1) %in% outliers_ACC_rna$outliers)]

expressionMatrix_ACC_rna <- read.csv('LUAD.csv', header = TRUE) # LOAD PROCESSED DATA
rownames(expressionMatrix_ACC_rna) <- expressionMatrix_ACC_rna[,1]
expressionMatrix_ACC_rna[,1] <- NULL
labels_ACC_rna <- expressionMatrix_ACC_rna[,dim(expressionMatrix_ACC_rna)[2]]
expressionMatrix_ACC_rna <- expressionMatrix_ACC_rna[,-dim(expressionMatrix_ACC_rna)[2]]
expressionMatrix_ACC_rna <- batchEffectRemoval(as.matrix(t(expressionMatrix_ACC_rna)), labels_ACC_rna, method = "sva") #BATCH EFFECT REMOVAL

#SCC
#countsInformation_SCC <- countsToMatrix('LUSC_pacientes.csv', extension = "")
#countsMatrix_SCC <- countsInformation_SCC$countsMatrix
#labels_SCC_rna <- countsInformation_SCC$labels 
#for(i in 1:length(labels_SCC_rna)){
  if (labels_SCC_rna[i] == 'Escamosas'){
    labels_SCC_rna[i] = 'SCC'
  }
  else {
    labels_SCC_rna[i] = 'Control'
  }
}
#myAnnotation_SCC <- getGenesAnnotation(rownames(countsMatrix_SCC))
#expressionMatrix_SCC_rna <- calculateGeneExpressionValues(countsMatrix_SCC,myAnnotation_SCC,genesNames = TRUE)
#outliers_SCC_rna <- RNAseqQA(expressionMatrix_SCC_rna, toRemoval = TRUE, toPNG = FALSE, toPDF = TRUE, outdir = 'D:/documentos/articulo')
#expressionMatrix_SCC_rna <- expressionMatrix_SCC_rna[,-which(colnames(expressionMatrix_SCC_rna) %in% outliers_SCC_rna$outliers)]
#labels_SCC_rna <- labels_SCC_rna[-which(colnames(expressionMatrix_SCC_rna) %in% outliers_SCC_rna$outliers)]

expressionMatrix_SCC_rna <- read.csv('LUSC.csv', header = TRUE) # LOAD PROCESSED DATA
rownames(expressionMatrix_SCC_rna) <- expressionMatrix_SCC_rna[,1]
expressionMatrix_SCC_rna[,1] <- NULL
labels_SCC_rna <- expressionMatrix_SCC_rna[,dim(expressionMatrix_SCC_rna)[2]]
expressionMatrix_SCC_rna <- expressionMatrix_SCC_rna[,-dim(expressionMatrix_SCC_rna)[2]]
expressionMatrix_SCC_rna <- batchEffectRemoval(as.matrix(t(expressionMatrix_SCC_rna)), labels_SCC_rna, method = "sva") #BATCH EFFECT REMOVAL

# ACC AND SCC

# MICROARRAY
expressionMatrix_tri_mic_1 <- mic_exp
labels_tri_mic_1 <- labels_tri_mic_1


# RNA-Seq
expressionMatrix_tri_rna <- read.csv('tri.csv', header = TRUE) # LOAD PROCESSED DATA
rownames(expressionMatrix_tri_rna) <- expressionMatrix_tri_rna[,1]
expressionMatrix_tri_rna[,1] <- NULL
labels_tri_rna <- expressionMatrix_tri_rna[,dim(expressionMatrix_tri_rna)[2]]
expressionMatrix_tri_rna <- expressionMatrix_tri_rna[,-dim(expressionMatrix_tri_rna)[2]]
expressionMatrix_tri_rna <- batchEffectRemoval(as.matrix(t(expressionMatrix_tri_rna)), labels_tri_rna, method = "sva") #BATCH EFFECT REMOVAL


#                                    SCC VS Control (microarray gene signature assess over RNA-SEQ)

folds <- 10
expressionMatrix_SCC_mic_huella_1 <- expressionMatrix_SCC_mic_1[,which(colnames(expressionMatrix_SCC_mic_1) %in% rownames(expressionMatrix_SCC_rna))]
cvIndex_SCC_mic_1 <- createDataPartition(labels_SCC_mic_1, p = .80, 
                                         list = FALSE, 
                                         times = folds) 
cvResults_SCC_mic_1 <- list()
cvDEGs_SCC_mic_1 <- list ()
for (i in seq_len(folds)){
  cvResults_SCC_mic_1[[i]] <- DEGsExtraction(t(expressionMatrix_SCC_mic_huella_1[cvIndex_SCC_mic_1[,i],]), as.factor(labels_SCC_mic_1[cvIndex_SCC_mic_1[,i]]), lfc=3, pvalue = 0.01, number = Inf, CV=FALSE)
  cvDEGs_SCC_mic_1[[i]] <- rownames(cvResults_SCC_mic_1[[i]]$DEG_Results$DEGs_Table)
}
huella_limma_SCC_mic_1 <- Reduce(f='intersect', cvDEGs_SCC_mic_1)

huella_limma_mrmr_SCC_mic_1 <- featureSelection(expressionMatrix_SCC_mic_1, labels_SCC_mic_1, huella_limma_SCC_mic_1, mode = 'mrmr')

Index_SCC_rna <- createDataPartition(labels_SCC_rna, p = .80, 
                                     list = FALSE, 
                                     times = 1)

train_labels_SCC_rna <- labels_SCC_rna[Index_SCC_rna]
test_labels_SCC_rna <- labels_SCC_rna[-Index_SCC_rna]
train_SCC_rna <- expressionMatrix_SCC_rna[,Index_SCC_rna]
test_SCC_rna <- expressionMatrix_SCC_rna[,-Index_SCC_rna]

knn_train_SCC_3A_1 <- knn_trn(t(train_SCC_rna), as.factor(train_labels_SCC_rna), names(huella_limma_mrmr_SCC_mic_1),5) 
knn_SCC_3A_1 <- knn_test(t(train_SCC_rna), as.factor(train_labels_SCC_rna), t(test_SCC_rna), as.factor(test_labels_SCC_rna), names(huella_limma_mrmr_SCC_mic_1), bestK =  knn_train_SCC_3A_1$bestK)


# PLOTS 

# TRAINING
plot(knn_train_SCC_3A_1$accuracyInfo$meanAccuracy[1:15], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.68,1.0), panel.first = grid(col='gray45'), cex.axis=1.3, cex.lab=1.3)
lines(knn_train_SCC_3A_1$sensitivityInfo$meanSensitivity[1:15], col='blue', lwd=2, lty=2)
lines(knn_train_SCC_3A_1$specificityInfo$meanSpecificity[1:15], col='#FF8B00', lwd=2, lty=4)
lines(knn_train_SCC_3A_1$F1Info$meanF1[1:15], col='red', lwd=2, lty=4)
legend(x=11.53 ,y =0.756, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

# CONFUSION MATRIX
dataPlot(knn_SCC_3A_1$cfMats[[4]]$table,labels  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)




#                                                 SCC VS CONTROL (RNA-SEQ gene signature assess over microarray)


expressionMatrix_SCC_rna_huella_1 <- expressionMatrix_SCC_rna[which(rownames(expressionMatrix_SCC_rna)%in% colnames(expressionMatrix_SCC_mic_1)),]
cvIndex_SCC_rna <- createDataPartition(labels_SCC_rna, p = .80, 
                                       list = FALSE, 
                                       times = folds) 

cvResults_SCC_rna_1 <- list()
cvDEGs_SCC_rna_1 <- list ()
for (i in seq_len(folds)){
  cvResults_SCC_rna_1[[i]] <- DEGsExtraction(expressionMatrix_SCC_rna_huella_1[,cvIndex_SCC_rna[,i]], as.factor(labels_SCC_rna[cvIndex_SCC_rna[,i]]), lfc=3.5, pvalue = 0.01, number = Inf, CV = FALSE)
  cvDEGs_SCC_rna_1[[i]] <- rownames(cvResults_SCC_rna_1[[i]]$DEG_Results$DEGs_Table)
}
huella_limma_SCC_rna_1 <- Reduce(f='intersect', cvDEGs_SCC_rna_1)
huella_limma_mrmr_SCC_rna_1 <- featureSelection(t(expressionMatrix_SCC_rna), labels_SCC_rna, huella_limma_SCC_rna_1, mode = 'mrmr')

Index_SCC_mic_1 <- createDataPartition(labels_SCC_mic_1, p = .80, 
                                       list = FALSE, 
                                       times = 1)
train_labels_SCC_mic_1 <- labels_SCC_mic_1[Index_SCC_mic_1]
test_labels_SCC_mic_1 <- labels_SCC_mic_1[-Index_SCC_mic_1]
train_SCC_mic_1 <- expressionMatrix_SCC_mic_1[Index_SCC_mic_1,]
test_SCC_mic_1 <- expressionMatrix_SCC_mic_1[-Index_SCC_mic_1,]

knn_train_SCC_3B_1 <- knn_trn(train_SCC_mic_1, as.factor(train_labels_SCC_mic_1), names(huella_limma_mrmr_SCC_rna_1),5) 
knn_SCC_3B_1 <- knn_test(train_SCC_mic_1, as.factor(train_labels_SCC_mic_1), test_SCC_mic_1, as.factor(test_labels_SCC_mic_1), names(huella_limma_mrmr_SCC_rna_1), bestK =  knn_train_SCC_3B_1$bestK)

# CONFUSION MATRIX
dataPlot(knn_SCC_3B_1$cfMats[[4]]$table,labels  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)



#                                      ACC vs SCC vs Control (Microarray gene signature assess over RNA-Seq) 


expressionMatrix_tri_mic_1_huella <- expressionMatrix_tri_mic_1[,which(colnames(expressionMatrix_tri_mic_1)%in% rownames(expressionMatrix_tri_rna))]
cvIndex_tri_mic_1 <- createDataPartition(labels_tri_mic_1, p = .80, 
                                         list = FALSE, 
                                         times = folds) 
cvResults_tri_mic_1 <- list()
cvDEGs_tri_mic_1 <- list ()
for (i in seq_len(folds)){
  cvResults_tri_mic_1[[i]] <- DEGsExtraction(t(expressionMatrix_tri_mic_1_huella[cvIndex_tri_mic_1[,i],]), as.factor(labels_tri_mic_1[cvIndex_tri_mic_1[,i]]), lfc=1, cov = 2, pvalue = 0.01, number = Inf, CV = FALSE)
  cvDEGs_tri_mic_1[[i]] <- rownames(cvResults_tri_mic_1[[i]]$DEG_Results$MulticlassLFC)
}
funcion <- function(x){
  h <-vector()
  i <- 1
  num <-0
  while (i <= length(x)){
    for (j in 1:length(x[[i]])){
      if (!(x[[i]][j] %in% h)){
        for (k in 1:length(x)){
          if (x[[i]][j] %in% x[[k]]){
            
            num = num +1
          }
          
        }
        if (num >=8){
          h <- c(h, x[[i]][j])
          
        }
        num=0
        
      }
    }
    i = i+1
    
  }
  h
}
h1 <- funcion(cvDEGs_tri_mic_1)
huella_limma_mrmr_tri_mic_1a <- featureSelection(expressionMatrix_tri_mic_1, labels_tri_mic_1, h1, mode = 'mrmr')

Index_tri_rna <- createDataPartition(labels_tri_rna, p = .80, 
                                     list = FALSE, 
                                     times = 1)

train_labels_tri_rna <- labels_tri_rna[Index_tri_rna]
test_labels_tri_rna <- labels_tri_rna[-Index_tri_rna]
train_tri_rna <- expressionMatrix_tri_rna[,Index_tri_rna]
test_tri_rna <- expressionMatrix_tri_rna[,-Index_tri_rna]

knn_train_tri_4A_1a <- knn_trn(t(train_tri_rna), as.factor(train_labels_tri_rna), names(huella_limma_mrmr_tri_mic_1a),5) # k=17
knn_tri_4A_1a <- knn_test(t(train_tri_rna), as.factor(train_labels_tri_rna), t(test_tri_rna), as.factor(test_labels_tri_rna), names(huella_limma_mrmr_tri_mic_1a), bestK =  knn_train_tri_4A_1a$bestK)

plot(knn_train_tri_4A_1a$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.15,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
lines(knn_train_tri_4A_1a$sensitivityInfo$meanSensitivity, col='blue', lwd=2, lty=2)
lines(knn_train_tri_4A_1a$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2, lty=4)
lines(knn_train_tri_4A_1a$F1Info$meanF1, col='red', lwd=2, lty=4)
legend(x=7.77 ,y =0.352, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

dataPlot(knn_tri_4A_1a$cfMats[[8]]$table,labels  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)


#                                       ACC vs SCC vs Control (RNA-Seq gene signature assess over Microarray)


expressionMatrix_tri_rna_huella_1 <- expressionMatrix_tri_rna[which(rownames(expressionMatrix_tri_rna)%in% colnames(expressionMatrix_tri_mic_1)),]
cvIndex_tri_rna <- createDataPartition(labels_tri_rna, p = .80, 
                                       list = FALSE, 
                                       times = folds) 
cvResults_tri_rna_1 <- list()
cvDEGs_tri_rna_1 <- list ()
for (i in seq_len(folds)){
  cvResults_tri_rna_1[[i]] <- DEGsExtraction(expressionMatrix_tri_rna_huella_1[,cvIndex_tri_rna[,i]], as.factor(labels_tri_rna[cvIndex_tri_rna[,i]]), lfc=2, cov = 2, pvalue = 0.01, number = Inf, CV = FALSE)
  cvDEGs_tri_rna_1[[i]] <- rownames(cvResults_tri_rna_1[[i]]$DEG_Results$MulticlassLFC)
}
huella_limma_tri_rna_1 <- Reduce(f='intersect', cvDEGs_tri_rna_1)
huella_limma_mrmr_tri_rna_1 <- featureSelection(t(expressionMatrix_tri_rna), labels_tri_rna, huella_limma_tri_rna_1, mode = 'mrmr')
Index_tri_mic_1 <- createDataPartition(labels_tri_mic_1, p = .80, 
                                       list = FALSE, 
                                       times = 1)

train_labels_tri_mic_1<- labels_tri_mic_1[Index_tri_mic_1]
test_labels_tri_mic_1 <- labels_tri_mic_1[-Index_tri_mic_1]
train_tri_mic_1 <- expressionMatrix_tri_mic_1[Index_tri_mic_1,]
test_tri_mic_1 <- expressionMatrix_tri_mic_1[-Index_tri_mic_1,]

knn_train_tri_4B_1 <- knn_trn(train_tri_mic_1, as.factor(train_labels_tri_mic_1), names(huella_limma_mrmr_tri_rna_1),5) 
knn_tri_4B_1 <- knn_test(train_tri_mic_1, as.factor(train_labels_tri_mic_1), test_tri_mic_1, as.factor(test_labels_tri_mic_1), names(huella_limma_mrmr_tri_rna_1), bestK =  knn_train_tri_4B_1$bestK)

dataPlot(knn_tri_4B_1$cfMats[[8]]$table,labels  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)


# BOXPLOTS

genes <-  c('SCEL','ZNF91','CMTM6','HNRNPC')
box_mic_scc <- as.data.frame(expressionMatrix_SCC_mic_huella_1[,which(colnames(expressionMatrix_SCC_mic_huella_1) %in% genes)])
box_mic_scc1 <- box_mic_scc %>% select(SCEL,CMTM6,ZNF91,HNRNPC)
box_rna_scc <- as.data.frame(t(expressionMatrix_SCC_rna_huella_1[which(rownames(expressionMatrix_SCC_rna_huella_1) %in% genes),]))
box_rna_scc1 <- box_rna_scc %>% select(SCEL,CMTM6,ZNF91,HNRNPC)
box_matrix <- rbind(box_mic_scc1, box_rna_scc1)
labels_SCC <- labels_SCC_mic_1
for (i in 1:length(labels_SCC)){
  if (labels_SCC[i]=='Control'){
    labels_SCC[i] <- 'Control_mic'
  } else {
    labels_SCC[i] <- 'SCC_mic'
  } 
}
labels_SCC_1 <- labels_SCC_rna
for (i in 1:length(labels_SCC_1)){
  if (labels_SCC_1[i]=='Control'){
    labels_SCC_1[i] <- 'Control_rna'
  } else {
    labels_SCC_1[i] <- 'SCC_rna'
  } 
}
labels_box_SCC <- c(labels_SCC,labels_SCC_1)
labels_box_SCC_a <- labels_box_SCC[-c(30,482,424)]
labels_box_SCC_a <- c(labels_box_SCC_a[1],labels_box_SCC[c(30,482,424)], labels_box_SCC_a[2:length(labels_box_SCC_a)])
box_matrix_a <- box_matrix[-c(30,482,424),]
box_matrix_a <- rbind(box_matrix_a[1,],box_matrix[c(30,482,424),], box_matrix_a[2:dim(box_matrix_a)[1],])
dataPlot(t(as.matrix(box_matrix_a)),labels_box_SCC_a,mode = "genesBoxplot",toPNG = FALSE, colours = c("darkred", "forestgreen", "darkorchid3", "dodgerblue3"))
#triclase
genes_tri <- c('ZNF91','ADAM23','TSPAN3','CADM1','DAB2','TMED5','SORT1','PERP')

labels_tri <- labels_tri_mic_1
for (i in 1:length(labels_tri)){
  if (labels_tri[i]=='Control'){
    labels_tri[i] <- 'Control_mic'
  } else if (labels_tri[i]=='SCC') {
    labels_tri[i] <- 'SCC_mic'
  } else {
    labels_tri[i] <- 'ACC_mic'
  } 
}

labels_tri_1 <- labels_tri_rna
for (i in 1:length(labels_tri_1)){
  if (labels_tri_1[i]=='Control'){
    labels_tri_1[i] <- 'Control_rna'
  } else if (labels_tri_1[i]=='SCC') {
    labels_tri_1[i] <- 'SCC_rna'
  } else {
    labels_tri_1[i] <- 'ACC_rna'
  } 
}

labels_box_tri <- c(labels_tri,labels_tri_1)

box_mic_tri <- as.data.frame(expressionMatrix_tri_mic_1_huella[,which(colnames(expressionMatrix_tri_mic_1_huella) %in% genes_tri)])
box_mic_tri1 <- box_mic_tri %>% select(ZNF91,ADAM23,TSPAN3,CADM1,DAB2,TMED5,SORT1,PERP)

box_rna_tri <- as.data.frame(t(expressionMatrix_tri_rna_huella_1[which(rownames(expressionMatrix_tri_rna_huella_1) %in% genes_tri),]))
box_rna_tri1 <- box_rna_tri %>% select(ZNF91,ADAM23,TSPAN3,CADM1,DAB2,TMED5,SORT1,PERP)

box_matrix_tri <- rbind(box_mic_tri1, box_rna_tri1)

labels_box_tri_a <- labels_box_tri[-c(815,849,896)]
labels_box_tri_a <- c(labels_box_tri_a[1:814],labels_box_tri[c(896,815,849)], labels_box_tri_a[815:length(labels_box_tri_a)])
box_matrix_tri_a <- box_matrix_tri[-c(896,815,849),]
box_matrix_tri_a <- rbind(box_matrix_tri_a[1:814,],box_matrix_tri[c(896,815,849),], box_matrix_tri_a[815:dim(box_matrix_tri_a)[1],])
dataPlot(t(as.matrix(box_matrix_tri_a)),labels_box_tri_a,mode = "genesBoxplot",toPNG = FALSE, colours = c("darkred", "forestgreen", "darkorchid3", "dodgerblue3"))



# AUC

# triclasS 

response <- as.factor(test_labels_tri_rna)
aucs <- rep(NA, length(levels(response))) # store AUCs
legendLabels <- as.character()
colours <- c('red','blue','green')

par(oma = c(5, 1, 0, 1))
plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab="Sensitivity",
     xlab="1 - Specificity",
     bty='n',
     cex.lab=1.3,
     cex.axis=1.3)

for (i in seq_along(levels(response))) {
  cur.class <- levels(response)[i]
  binaryTraining.labels <- as.factor(train_labels_tri_rna == cur.class)
  binaryTest.labels <- as.factor(test_labels_tri_rna == cur.class)
  
  binary_knn_cv_mrmr_results <- knn_trn(t(train_tri_rna), binaryTraining.labels, names(huella_limma_mrmr_tri_mic_1a))
  
  binary_knn_test_mrmr_results <- knn_test(t(train_tri_rna), binaryTraining.labels, t(test_tri_rna), binaryTest.labels, names(huella_limma_mrmr_tri_mic_1a), bestK = binary_knn_cv_mrmr_results$bestK)
  
  score <- binary_knn_test_mrmr_results$predictions[[8]]
  score <- as.vector(score)
  score[score=='FALSE'] <- 0
  score[score=='TRUE'] <- 1
  binaryTest.labels <- as.vector(binaryTest.labels)
  binaryTest.labels[binaryTest.labels=='FALSE'] <- 0
  binaryTest.labels[binaryTest.labels=='TRUE'] <- 1
  pred <- prediction(as.numeric(score), as.numeric(binaryTest.labels))
  perf <- performance(pred, "tpr", "fpr")
  roc.x <- unlist(perf@x.values)
  roc.y <- unlist(perf@y.values)
  lines(roc.y ~ roc.x, col = colours[i], lwd = 2)
  # store AUC
  auc <- performance(pred, "auc")
  auc <- unlist(slot(auc, "y.values"))
  aucs[i] <- auc
  legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
}

print(paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))

lines(x=c(0,1), c(0,1))
par(fig = c(0, 1, 0, 1), oma = c(0.6, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)
