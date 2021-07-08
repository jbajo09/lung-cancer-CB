library(KnowSeq)
library(CORElearn)
library(dplyr)
library(caret)
library(car)
library(R.utils)
require(ROCR)
require(annotate)

set.seed(30)

# DATA READ AND PREPARATION

# MICROARRAY
mic_exp <- read.csv('mic_exp1.csv', header = TRUE)
mic_exp[,1] <- NULL
samples_mic <- read.csv('mic_samples.csv')
samples_mic[,1] <- NULL
expressionMatrix_tri_mic_1 <- as.matrix(mic_exp)
labels_tri_mic_1 <- samples_mic[,2]
outliers <- RNAseqQA(t(expressionMatrix_tri_mic_1), toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE, outdir = '.') 

#RNA-SEQ
expressionMatrix_tri_rna <- read.csv('tri.csv', header = TRUE) # LOAD PROCESSED DATA
rownames(expressionMatrix_tri_rna) <- expressionMatrix_tri_rna[,1]
expressionMatrix_tri_rna[,1] <- NULL
labels_tri_rna <- expressionMatrix_tri_rna[,dim(expressionMatrix_tri_rna)[2]]
expressionMatrix_tri_rna <- expressionMatrix_tri_rna[,-dim(expressionMatrix_tri_rna)[2]]
expressionMatrix_tri_rna <- expressionMatrix_tri_rna[,order(colnames(expressionMatrix_tri_rna))]
expressionMatrix_tri_rna <- expressionMatrix_tri_rna[,which(colnames(expressionMatrix_tri_rna) %in% colnames(expressionMatrix_tri_mic_1))]
expressionMatrix_tri_mic_1 <- expressionMatrix_tri_mic_1[,which(colnames(expressionMatrix_tri_mic_1) %in% colnames(expressionMatrix_tri_rna))]


#batch effect treatment
matrix <- rbind(expressionMatrix_tri_mic_1, expressionMatrix_tri_rna)
labels <- c(labels_tri_mic_1,labels_tri_rna)
matrix <- batchEffectRemoval(as.matrix(t(matrix)), labels, method = "sva") #BATCH EFFECT REMOVAL
expressionMatrix_tri_mic_1_2 <- matrix[,1:814]
expressionMatrix_tri_rna_2 <- matrix[,815:dim(matrix)[2]]

# SCC vs control matrix
#rna
LUSC_samples <- read.csv('LUSC_pacientes.csv')
LUSC_samples <- LUSC_samples$Run
labels_SCC_rna <- labels_tri_rna[which(colnames(expressionMatrix_tri_rna_2)%in% LUSC_samples)]
expressionMatrix_SCC_rna <- expressionMatrix_tri_rna_2[,which(colnames(expressionMatrix_tri_rna_2)%in% LUSC_samples)]
#mic
labels_SCC_mic <- labels_tri_mic_1[which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'SCC')]
expressionMatrix_SCC_mic <- expressionMatrix_tri_mic_1_2[,which(labels_tri_mic_1 == 'Control' | labels_tri_mic_1 == 'SCC')]



##                                  ACC VS SCC VS CONTROL

#                                   HUELLA MICROARRAY
# 
set.seed(2)
DEGs_CV_2 <- DEGsExtraction(expressionMatrix_tri_mic_1_2, as.factor(labels_tri_mic_1), lfc=1, cov = 2, pvalue = 0.05, number = Inf, CV = TRUE, numFolds=10)
DEGs_2 <- DEGs_CV_2$Common_DEGs
huella_limma_mrmr_tri_mic_1a_2 <- featureSelection(t(expressionMatrix_tri_mic_1_2), labels_tri_mic_1, DEGs_2, mode = 'mrmr')

set.seed(333)
Index_tri_rna1 <- createDataPartition(labels_tri_rna, p = .80, 
                                     list = FALSE, 
                                     times = 1)

train_labels_tri_rna_2 <- labels_tri_rna[Index_tri_rna]
test_labels_tri_rna_2 <- labels_tri_rna[-Index_tri_rna]
train_tri_rna_2 <- expressionMatrix_tri_rna_2[,Index_tri_rna]
test_tri_rna_2 <- expressionMatrix_tri_rna_2[,-Index_tri_rna]
knn_train_tri_4A_1a_2 <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), names(huella_limma_mrmr_tri_mic_1a_2),5) # k=7
knn_tri_4A_1a_2 <- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2), names(huella_limma_mrmr_tri_mic_1a_2), bestK =  knn_train_tri_4A_1a_2$bestK)


plot(knn_train_tri_4A_1a_2$accuracyInfo$meanAccuracy[1:20], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.79,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
lines(knn_train_tri_4A_1a_2$sensitivityInfo$meanSensitivity[1:20], col='blue', lwd=2, lty=2)
lines(knn_train_tri_4A_1a_2$specificityInfo$meanSpecificity[1:20], col='#FF8B00', lwd=2, lty=4)
lines(knn_train_tri_4A_1a_2$F1Info$meanF1[1:20], col='red', lwd=2, lty=4)
legend(x=15.9 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

dataPlot(knn_tri_4A_1a_2$cfMats[[8]]$table,labels_tri_rna  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)


response <- as.factor(test_labels_tri_rna_2)
aucs <- rep(NA, length(levels(response))) # store AUCs
legendLabels <- as.character()
colours <- c('red','blue','green')


plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab="Sensitivity",
     xlab="1 - Specificity",
     bty='n',
     cex.lab=1.3,
     cex.axis=1.3)
for (i in seq_along(levels(response))) {
  cur.class <- levels(response)[i]
  binaryTraining.labels <- as.factor(train_labels_tri_rna_2 == cur.class)
  binaryTest.labels <- as.factor(test_labels_tri_rna_2 == cur.class)
  
  tri_knn_cv_mrmr_results <- knn_trn(t(train_tri_rna_2), binaryTraining.labels, names(huella_limma_mrmr_tri_mic_1a_2))
  
  tri_knn_test_mrmr_results <- knn_test(t(train_tri_rna_2), binaryTraining.labels, t(test_tri_rna_2), binaryTest.labels, names(huella_limma_mrmr_tri_mic_1a_2), bestK = tri_knn_cv_mrmr_results$bestK)
  
  score <- tri_knn_test_mrmr_results$predictions[[8]]
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

print(paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2))) #0.97

lines(x=c(0,1), c(0,1))
legend(x=0.48 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)


#triclase
genes_tri <- c('KRT15','TFRC','S100A2','COL17A1','GPX2','FERMT1','SPRR1B','ZNF750')

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

box_mic_tri <- as.data.frame(t(expressionMatrix_tri_mic_1_2[which(rownames(expressionMatrix_tri_mic_1_2) %in% genes_tri),]))
box_mic_tri <- box_mic_tri %>% select(KRT15,TFRC,S100A2,COL17A1,GPX2,FERMT1,SPRR1B,ZNF750)

box_rna_tri <- as.data.frame(t(expressionMatrix_tri_rna_2[which(rownames(expressionMatrix_tri_rna_2) %in% genes_tri),]))
box_rna_tri <- box_rna_tri %>% select(KRT15,TFRC,S100A2,COL17A1,GPX2,FERMT1,SPRR1B,ZNF750)

box_matrix_tri <- rbind(box_mic_tri, box_rna_tri)

labels_box_tri_a <- labels_box_tri[-c(815,849,896)]
labels_box_tri_a <- c(labels_box_tri_a[1:814],labels_box_tri[c(896,815,849)], labels_box_tri_a[815:length(labels_box_tri_a)])
box_matrix_tri_a <- box_matrix_tri[-c(896,815,849),]
box_matrix_tri_a <- rbind(box_matrix_tri_a[1:814,],box_matrix_tri[c(896,815,849),], box_matrix_tri_a[815:dim(box_matrix_tri_a)[1],])
dataPlot(t(as.matrix(box_matrix_tri_a)),labels_box_tri_a,mode = "genesBoxplot",toPNG = FALSE, colours = c("tomato", "chartreuse", "cyan", "red",'green4','blue'),colnumber = 3)


#                                   HUELLA RNA-SEQ
#
DEGs_CV_rna_2 <- DEGsExtraction(expressionMatrix_tri_rna_2, as.factor(labels_tri_rna), lfc=3, cov = 2, pvalue = 0.05, number = Inf, CV = TRUE, numFolds=10)
DEGs_rna_2 <- DEGs_CV_rna_2$Common_DEGs
huella_limma_tri_rna_1_2 <- featureSelection(t(expressionMatrix_tri_rna_2), labels_tri_rna, DEGs_rna_2, mode = 'mrmr')
Index_tri_mic_1 <- createDataPartition(labels_tri_mic_1, p = .80, 
                                       list = FALSE, 
                                       times = 1)

train_labels_tri_mic_1_2<- labels_tri_mic_1[Index_tri_mic_1]
test_labels_tri_mic_1_2 <- labels_tri_mic_1[-Index_tri_mic_1]
train_tri_mic_1_2 <- expressionMatrix_tri_mic_1_2[,Index_tri_mic_1]
test_tri_mic_1_2 <- expressionMatrix_tri_mic_1_2[,-Index_tri_mic_1]
knn_train_tri_4B_1_2 <- knn_trn(t(train_tri_mic_1_2), as.factor(train_labels_tri_mic_1_2), names(huella_limma_tri_rna_1_2)[1:20],5) # k=9
knn_tri_4B_1_2 <- knn_test(t(train_tri_mic_1_2), as.factor(train_labels_tri_mic_1_2), t(test_tri_mic_1_2), as.factor(test_labels_tri_mic_1_2), names(huella_limma_tri_rna_1_2)[1:20], bestK =  knn_train_tri_4B_1_2$bestK)

dataPlot(knn_tri_4B_1_2$cfMats[[8]]$table,labels_tri_mic_1  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)

plot(knn_train_tri_4B_1_2$accuracyInfo$meanAccuracy[], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.4,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
lines(knn_train_tri_4B_1_2$sensitivityInfo$meanSensitivity[], col='blue', lwd=2, lty=2)
lines(knn_train_tri_4B_1_2$specificityInfo$meanSpecificity[], col='#FF8B00', lwd=2, lty=4)
lines(knn_train_tri_4B_1_2$F1Info$meanF1[], col='red', lwd=2, lty=4)
legend(x=15.9 ,y =0.544, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

response <- as.factor(test_labels_tri_mic_1_2)
aucs <- rep(NA, length(levels(response))) # store AUCs
legendLabels <- as.character()
colours <- c('red','blue','green')


plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab="Sensitivity",
     xlab="1 - Specificity",
     bty='n',
     cex.lab=1.3,
     cex.axis=1.3)

for (i in seq_along(levels(response))) {
  cur.class <- levels(response)[i]
  binaryTraining.labels <- as.factor(train_labels_tri_mic_1_2 == cur.class)
  binaryTest.labels <- as.factor(test_labels_tri_mic_1_2 == cur.class)
  
  tri_knn_cv_mrmr_results <- knn_trn(t(train_tri_mic_1_2), binaryTraining.labels, names(huella_limma_tri_rna_1_2))
  
  tri_knn_test_mrmr_results <- knn_test(t(train_tri_mic_1_2), binaryTraining.labels, t(test_tri_mic_1_2), binaryTest.labels, names(huella_limma_tri_rna_1_2), bestK = tri_knn_cv_mrmr_results$bestK)
  
  score <- tri_knn_test_mrmr_results$predictions[[9]]
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
legend(x=0.48 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)



#                                        SCC VS CONTROL
#                                     HUELLA MICROARRAY
DEGs_CV_SCC_mic <- DEGsExtraction(expressionMatrix_SCC_mic, as.factor(labels_SCC_mic), lfc=4, pvalue = 0.01, number = Inf, CV = TRUE, numFolds=10)
DEGs_SCC_mic <- DEGs_CV_SCC_mic$Common_DEGs
huella_limma_mrmr_SCC_mic <- featureSelection(t(expressionMatrix_SCC_mic), labels_SCC_mic, DEGs_SCC_mic, mode = 'mrmr')
Index_SCC_rna <- createDataPartition(labels_SCC_rna, p = .80, 
                                     list = FALSE, 
                                     times = 1)

train_labels_SCC_rna <- labels_SCC_rna[Index_SCC_rna]
test_labels_SCC_rna <- labels_SCC_rna[-Index_SCC_rna]
train_SCC_rna <- expressionMatrix_SCC_rna[,Index_SCC_rna]
test_SCC_rna <- expressionMatrix_SCC_rna[,-Index_SCC_rna]

knn_train_SCC_rna <- knn_trn(t(train_SCC_rna), as.factor(train_labels_SCC_rna), names(huella_limma_mrmr_SCC_mic),5) # k=7
knn_SCC_rna <- knn_test(t(train_SCC_rna), as.factor(train_labels_SCC_rna), t(test_SCC_rna), as.factor(test_labels_SCC_rna), names(huella_limma_mrmr_SCC_mic), bestK =  knn_train_SCC_rna$bestK)


plot(knn_train_SCC_rna$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.79,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
lines(knn_train_SCC_rna$sensitivityInfo$meanSensitivity, col='blue', lwd=2, lty=2)
lines(knn_train_SCC_rna$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2, lty=4)
lines(knn_train_SCC_rna$F1Info$meanF1, col='red', lwd=2, lty=4)
legend(x=4.92 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

dataPlot(knn_SCC_rna$cfMats[[4]]$table,labels  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)



#BOXPLOT
genes <- c("COL1A1", "KRT15",  "SPRR1B", "KRT19")

box_mic_scc <- as.data.frame(t(expressionMatrix_SCC_mic[which(rownames(expressionMatrix_SCC_mic) %in% genes),]))
box_mic_scc1 <- box_mic_scc %>% select(COL1A1, KRT15,  SPRR1B, KRT19)
box_rna_scc <- as.data.frame(t(expressionMatrix_SCC_rna[which(rownames(expressionMatrix_SCC_rna) %in% genes),]))
box_rna_scc1 <- box_rna_scc %>% select(COL1A1, KRT15,  SPRR1B, KRT19)
box_matrix <- rbind(box_mic_scc1, box_rna_scc1)
labels_SCC <- labels_SCC_mic
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
dataPlot(t(as.matrix(box_matrix_a)),labels_box_SCC_a,mode = "genesBoxplot",toPNG = FALSE, colours = c("darkred", "forestgreen", "darkorchid3", "dodgerblue3"),colnumber = 2)

#                                   HUELLA RNA-SEQ
# NEW 2
DEGs_CV_SCC_rna <- DEGsExtraction(expressionMatrix_SCC_rna, as.factor(labels_SCC_rna), lfc=9, pvalue = 0.01, number = Inf, CV = TRUE, numFolds=10)
DEGs_SCC_rna <- DEGs_CV_SCC_rna$Common_DEGs
huella_limma_mrmr_SCC_rna <- featureSelection(t(expressionMatrix_SCC_rna), labels_SCC_rna, DEGs_SCC_rna, mode = 'mrmr')
Index_SCC_mic <- createDataPartition(labels_SCC_mic, p = .80, 
                                       list = FALSE, 
                                       times = 1)

train_labels_SCC_mic<- labels_SCC_mic[Index_SCC_mic]
test_labels_SCC_mic <- labels_SCC_mic[-Index_SCC_mic]
train_SCC_mic <- expressionMatrix_SCC_mic[,Index_SCC_mic]
test_SCC_mic <- expressionMatrix_SCC_mic[,-Index_SCC_mic]
knn_train_SCC_mic <- knn_trn(t(train_SCC_mic), as.factor(train_labels_SCC_mic), names(huella_limma_mrmr_SCC_rna),5) # k=5
knn_SCC_mic <- knn_test(t(train_SCC_mic), as.factor(train_labels_SCC_mic), t(test_SCC_mic), as.factor(test_labels_SCC_mic), names(huella_limma_mrmr_SCC_rna), bestK =  knn_train_SCC_mic$bestK)

plot(knn_train_SCC_mic$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.79,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
lines(knn_train_SCC_mic$sensitivityInfo$meanSensitivity, col='blue', lwd=2, lty=2)
lines(knn_train_SCC_mic$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2, lty=4)
lines(knn_train_SCC_mic$F1Info$meanF1, col='red', lwd=2, lty=4)
legend(x=19.92 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)

dataPlot(knn_SCC_mic$cfMats[[4]]$table,labels_SCC_mic  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)



###                                      COMPARING GENE SIGNATURES                                 ###



#                                             SCC VS ACC VS CONTROL                                   #

# 0.94
# 0.95


# GDC 
huella_GDC <- c('TP53','TTN','MUC16','CSMD3','RYR2','USH2A','LRP1B','ZFHX4','XIRP2','SYNE1','SPTA1','FLG')
knn_train_huella_GDC_tri <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_GDC[c(1,2,5,6,7,8,10,11,12)],5) 
knn_huella_GDC_tri <- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),  huella_GDC[c(1,2,5,6,7,8,10,11,12)], bestK =  knn_train_huella_GDC_tri$bestK)
#0.83
#0.84

# Analysis of gene expression profiles of lung cancer subtypes with machine learning algorithms
huella_FeiYuan  <- c('CSTA', 'TP63', 'SERPINB13', 'CLCA2', 'BICD2', 'PERP', 'FAT2', 'BNC1', 'ATP11B', 'FAM83B', 'KRT5', 'PARD6G', 'PKP1')
knn_train_huella_FeiYuan <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_FeiYuan[1:9],5) 
knn_huella_huella_FeiYuan <- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),  huella_FeiYuan[1:9], bestK =  knn_train_huella_FeiYuan$bestK)
#0.93
#0.92




# Adaptive Unsupervised Feature Learning for Gene Signature Identification in Non-Small-Cell Lung Cancer
huella_XIUCAI_YE <- c('KRT6A','SFTPA2','SFTPA1','KRT14','PI3','GAPDH', 'S100A2', 'KRT16', 'SPRR1B', 'SPRR1A','CEACAM6', 'KRT6B' ,'PERP', 'S100A7','AKR1B10','PKP1','CSTA')
knn_train_huella_XIUCAI_YE <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_XIUCAI_YE[4:17],5) 
knn_huella_huella_XIUCAI_YE <- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),  huella_XIUCAI_YE[4:17], bestK =  knn_train_huella_XIUCAI_YE$bestK)
#0.94
#0.94

# Comparison of lung cancer cell lines representing four histopathological subtypes with gene expression profiling using quantitative real-time PCR
huella_Takashi <- c('CDH1', 'FOXG1', 'IGSF3', 'ISL1', 'MALL', 'PLAU', 'RAB25', 'S100P', 'SLCO4A1', 'STMN1', 'TGM2')
knn_train_huella_Takashi <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_Takashi,5) 
knn_huella_huella_Takashi<- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),  huella_Takashi, bestK =  knn_train_huella_Takashi$bestK)
# 0.85
# 0.87


# Identification of expression signatures for non-small-cell lung carcinoma subtype classification
huella_Ran <- c('CALML3','CERS3','DSG3','KRT16','KRT6A','KRT6B','LINC01206','SPRR1B','BUB1B','NCAPH','TMPRSS11A','PERP','KRT31')
knn_train_huella_Ran <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_Ran[which(huella_Ran %in% rownames(expressionMatrix_tri_rna))],5) 
knn_huella_huella_Ran<- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),   huella_Ran[which(huella_Ran %in% rownames(expressionMatrix_tri_rna))], bestK =  knn_train_huella_Ran$bestK)
#0.93
#0.90

# An Expression Signature as an Aid to the Histologic Classification of Non-Small Cell Lung Cancer
huella_Luc <- c('SPINK1', 'SFTA2','NKX2-1','CAPN8','CLDN3','TESC','LGSN','HNF1B','GOLT1A','HPN','FMO5','ACSL5','CDH15','ALDH3B1','RORC','DPP4','PRR15L','STK32A',
                'KCNK5','ABCC6','SMPDL3B','FAM83B','DST','COL7A1','SOX15','PNCK','VSNL1','DSC3','KRT16','SERPINB5','DAPL1','KRT17','PKP1','DSG3',
                'KRT13','CLCA2','KRT6C','TP63','KRT6B','CALML3','KRT6A','KRT5')
knn_train_huella_Luc <- knn_trn(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), huella_Luc[which(huella_Luc %in% rownames(expressionMatrix_tri_rna))],5) 
knn_huella_huella_Luc<- knn_test(t(train_tri_rna_2), as.factor(train_labels_tri_rna_2), t(test_tri_rna_2), as.factor(test_labels_tri_rna_2),   huella_Luc[which(huella_Luc %in% rownames(expressionMatrix_tri_rna))], bestK =  knn_train_huella_Luc$bestK)

#0.95   (0.85)
#0.95   (0.84)














