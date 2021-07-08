# Before the deploy to shinyapps.io, run:
# library(BiocManager)
# options(repos = BiocManager::repositories())
# See here:
# https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/5


# Packages
library(shiny)
library(shinydashboard)
library(dplyr)
library(KnowSeq)
library(reshape2)
library(caret)
library(ggplot2)
library(ggalluvial)
library(DT)
library(waiter)
library(ROCR)

# File to slightly modify dataPlot function
#source("www/dataPlot.R")

# Define some spinners
spinner_abrir <- tagList(
  spin_folding_cube(),
  span(br(), h4("Loading application..."), style="color:white;")
)

spinner <- tagList(
  spin_chasing_dots(),
  span(br(), h4("Loading..."), style="color:white; display: inline-block;")
)

ui <- dashboardPage(title = "KnowSeq integration", # Title in web browser
  ## Theme
  skin = "black",
  ## Header
  dashboardHeader(title = span(
    "KnowSeq",
    style = "font-family: Lucida Console; font-weight: bold"
  )),
  ## Sidebar
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("file-alt")),
      menuItem("Data loading", tabName = "datos", icon = icon("database")),
      menuItem("Genes selection", tabName = "genes", icon = icon("dna")),
      menuItem("Model training", tabName = "entrenamiento", icon = icon("play")),
      menuItem("Model validation", tabName = "validation", icon = icon("check-circle")),
      menuItem("Code", tabName = "codigo", icon = icon("code"))
    )
  ),
  ## Body
  dashboardBody(
    use_waiter(),
    # Spinners to show on load, or when the application is busy
    #waiter_show_on_load(spinner_abrir, color = "#027368"),
    #waiter_on_busy(spinner, color = "#027368"),
    tabItems(
      # Tab 1
      tabItem(tabName = "intro",
              
              h1("About this web application"),
              tags$i("lung cancer biomarkeRs"), "is a web application that allows users with no previous knowledge of programming to replicate integration between Microarray and RNA-Seq carried out in Heterogeneous Gene Expression Cross-Evaluation of Robust Biomarkers
Using Machine Learning Techniques Applied to Lung Cancer ",
              br(), br(),
              
              "The ", tags$i("lung cancer biomarkeRs"), "application is part of the", 
              tags$a(
                " work entitled Heterogeneous Gene Expression Cross-Evaluation of Robust Biomarkers
Using Machine Learning Techniques Applied to Lung Cancer",
                href = "https://github.com/jbajo09/lung-cancer-CB.git",
                target="_blank"
              ),
              "It's developed in R-Shiny and the code is ",
              tags$a(
                "open source.",
                href = "https://github.com/jbajo09/lung-cancer-CB.git",
                target="_blank"
              ),
              
              h2("Abstract "),
              
              h3(tags$b("Background")),
              "Nowadays, gene expression analysis is one of the most promising pillars for
                understanding and uncovering the mechanisms underlying the development and spread of cancer . In
                this sense, Next Generation Sequencing technologies, such as RNA-Seq, are currently leading the
                market due to their precision and cost. Nevertheless, there is still an enormous amount of
                non-analyzed data obtained from older technologies, such as Microarray, which could still be useful to
                extract relevant knowledge.",
              h3(tags$b("Methods")),
              "Throughout this research, a complete machine learning methodology to cross-evaluate the compatibility between both RNA-Seq 
              and Microarray sequencing technologies is described and implemented. In order to show a real application of the designed pipeline, 
              a lung cancer case study is addressed by considering two detected subtypes: adenocarcinoma and squamous cell carcinoma. 
              Transcriptomic datasets considered for our study have been obtained from the public repositories NCBI/GEO, ArrayExpress and GDC-Portal. 
              From them, several gene experiments have been carried out with the aim of finding gene signatures for these lung cancer subtypes, 
              linked to both transcriptomic technologies. With these DEGs selected, intelligent predictive models capable of classifying new samples 
              belonging to these cancer subtypes were developed. ",
              h3(tags$b("Results")),
              "The predictive models built using one technology are capable of discerning samples from a different technology. The classification results 
              are evaluated in terms of accuracy, F1-score and ROC curves along with AUC. Finally, the biological information of the gene sets obtained 
              and their relationship with lung cancer is reviewed, encountering strong biological evidence linking them to the disease. ",
            
              h3(tags$b("Conclusions")),
              "Our method has the capability of finding strong gene signatures which are also independent to the transcriptomic technology used to develop the analysis. 
              In addition, our article highlighted the potential of using heterogeneous transcriptomic data to increase the amount of samples for the studies, increasing 
              the statistical significance of the results.",
              

              # Images
              br(), br(), br(),
              fluidRow(column(6, tags$img(src = "ugr.png", height = "100px")),
                       column(6, tags$img(src = "knowseq.png", height = "120px")))
              
      ),
      
      # Tab 2
      tabItem(tabName = "datos",

              # Left column
              fluidRow(column(6, 
                h1("Data loading"),
                fileInput(inputId = "file_labels_mic",
                          label = span("Select CSV file with Microarray labels (see ",
                                       tags$a(
                                         "here",
                                         href = "https://drive.google.com/file/d/1oldhrr6I9ES7cyqLl37RuL13yf63cvvS/view?usp=sharing",
                                         target="_blank"
                                         ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                fileInput(inputId = "file_DEGsMatrix_mic",
                          label = span("Select CSV file with Microarray expression matrix (see ",
                                       tags$a(
                                         "here",
                                         href = "https://drive.google.com/file/d/1-zsAU3uiqB-0qspmL303Pbv5NF5U1G0v/view?usp=sharing",
                                         target="_blank"
                                       ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                fileInput(inputId = "file_labels_rna",
                          label = span("Select CSV file with RNA-Seq labels (see ",
                                       tags$a(
                                         "here",
                                         href = "https://drive.google.com/file/d/1hrNfPchAgVO_5ZNLHGrzS-PZr8PChvTu/view?usp=sharing",
                                         target="_blank"
                                       ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                fileInput(inputId = "file_DEGsMatrix_rna",
                          label = span("Select CSV file with RNA-Seq expression matrix (see ",
                                       tags$a(
                                         "here",
                                         href = "https://drive.google.com/file/d/1lhlgr1wXhwniOwzuOYPTHa6dPCk6R2TA/view?usp=sharing",
                                         target="_blank"
                                       ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                
                actionButton(inputId = "boton_importar",
                             label = "Import file",
                             icon = icon("fas fa-file-import", lib = "font-awesome"),
                             width = "100%"),
                br(),
                
                conditionalPanel(condition = "input.boton_importar!=0",
                                 
                                 h2("Distribution of classes in Microarray data"),
                                 
                                 tableOutput("tabla1")),
                br(),
                 
                conditionalPanel(condition = "input.boton_importar!=0",
                                 
                                 h2("Distribution of classes in RNA-Seq data"),
                                 
                                 tableOutput("tabla2"))
 
              
      
                            
                                 
                                 
                                 ),
                 # Right column
                 column(6, br(), br(), br(),
                          conditionalPanel(condition = "input.boton_importar!=0",
                             h2("Train-test partition"),
                             
                             sliderInput("porcentaje_entrenamiento",
                                         label = "Train percentage (%)",
                                         value = 80, min = 5, max = 95, step = 5,
                                         width = "100%"
                             ),
                           ))
              )

      ),
      
      # Tab 3
      tabItem(tabName = "genes",
              h1("DEGs extraction 10-folds CV"),
              textInput(inputId = "LFC", label = "Select the threshold LFC", value = 1, width = "50%"),
              
              sliderInput(inputId = "COV", label = "Select the threshold COV", value = 2, min = 1, max = 3, step = 1, width = "50%"),
              
              textInput(inputId = "pvalue", label = "Select the threshold p-value", value = 0.05, width = "50%"),
              
              actionButton(inputId = "boton_genes",
                           label = "Extract DEGs",
                           icon = icon("dna", lib = "font-awesome"),
                           width = "50%"),
              br(),
              br(),
              conditionalPanel(condition = "input.boton_genes!=0",
                                 
                h3("Number of DEGs extracted"),
                br(),
                textOutput("numberDEGs"),
                br(),
                fluidRow(
                  column(4, h4(tags$b("  MRMR")), tableOutput("genes_mrmr")),
                
                )
              )
      ),
      
      # Tab 4
      tabItem(tabName = "entrenamiento",
              h1("Model training"),
              
              # Choose number of folds
              selectInput("number_folds",
                          label = "Number of folds",
                          choices = c(3, 5, 10),
                          selected = 5,
                          width = "50%"),
              
              # Train model button
              actionButton(inputId = "boton_model_training",
                           label 
                           = "Train model",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              
              br(),
              br(),
              # OPtimal k
              textOutput("optimal_knn"),

              br(),
              h2("Training CV plot"),
              br(),
              plotOutput("train", width = "90%", height = '600px'),
              br(),
              br(),
              h2('Gene expression Boxplot'),
              actionButton(inputId = "boton_boxplot",
                           label 
                           = "BoxPlots",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              textInput(inputId = "nboxplot", label = "Select the number of genes from mRMR ranking to get boxplot", value = 9, width = "50%"),
              plotOutput("boxplot", width = "90%", height = '600px')
              
              
      ),

      # Tab 5
      tabItem(tabName = "validation",
              h1("Model validation"),
              
              
              sliderInput(inputId = "numero_genes_validation", label = "Select the number of genes to use:",
                          value = 9, min = 1, max = 50, step = 1, width = "50%"),

              actionButton(inputId = "boton_model_validation",
                           label = "Validate model in test",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              
              br(),
              br(),
              
              plotOutput("results_validation",
                         width = "80%"),
              br(),
              br(),
              
              h2('Area Under the Curve'),
              actionButton(inputId = "boton_calculate_AUC",
                           label = "Calculate AUC",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              textOutput("meanAUC"),
              plotOutput("AUC",
                         width = "80%"),
              
      ),
      
      # Tab 6
      tabItem(tabName = "codigo",
              h1("Code"),
              tags$h4(
                "In ", tags$a(href = "https://drive.google.com/drive/folders/1hH-V2F99mkic5BIKV1Plh6kDMcTyRBiF?usp=sharing", "this repository"),
                "you can find the code of this web application.")
              )
      ) # Close tabs
  ) # Close dashboard body
) # Close dashboard page

# Extend size of accepted files (40MB instead of the 5MB - default)
options(shiny.maxRequestSize = 600*1024^2)

server <- function(input, output){

  values <- reactiveValues(ranking = NULL, DEGsMatrix_rna_batch = NULL, optimalSVM_train = NULL, optimalkNN_train = NULL)
  
  # Server of tab: Data loading ------
  
  observeEvent(input$boton_importar, {
    
    # If files are selected, they are imported
    # Read labels
    labels_mic <- as.vector(t(read.csv2(file = input$file_labels_mic$datapath)))
    labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
    
    
    # Table
    output$tabla1 <- renderTable({
        if(is.null(input$file_labels_mic$datapath)) return(NULL)
        
      # Message if file is correctly imported
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " File imported"),
        easyClose = TRUE,
        footer = NULL
      ))
    
      tabla_aux <- as.data.frame(table(labels_mic)) %>% rename(Label = labels_mic, Samples = Freq)
      return(tabla_aux)
    })
    
    output$tabla2 <- renderTable({
      if(is.null(input$file_labels_rna$datapath)) return(NULL)
      
      # Message if file is correctly imported
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " Files imported"),
        easyClose = TRUE,
        footer = NULL
      ))
      
      tabla_aux1 <- as.data.frame(table(labels_rna)) %>% rename(Label = labels_rna, Samples = Freq)
      return(tabla_aux1)
    })
  
  }) # Close import button
  
  
  # Server of tab: Genes selection ------
  
  w <- Waiter$new(html = tagList(spin_folding_cube(),
                                span(br(), br(), br(), h4("Running DEGs 10-folds CV..."),
                                     style="color:white;")))
  
  observeEvent(input$boton_genes, {
    
    w$show()
    
    # Read labels
    labels_mic <- as.vector(t(read.csv2(file = input$file_labels_mic$datapath)))
    #labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
    # Read DEGsMatrix
    DEGsMatrix_mic <- as.matrix(read.csv(file = input$file_DEGsMatrix_mic$datapath))
    #DEGsMatrix_rna <- as.matrix(read.csv(file = input$file_DEGsMatrix_rna$datapath))
    
    #DEGsMatrix_mic <- DEGsMatrix_mic[,which(colnames(DEGsMatrix_mic) %in% colnames(DEGsMatrix_rna))]
    #DEGsMatrix_rna <- DEGsMatrix_rna[,which(colnames(DEGsMatrix_rna) %in% colnames(DEGsMatrix_mic))]
    #DEGsMatrix_mic <- DEGsMatrix_mic[,order(colnames(DEGsMatrix_mic))]
    #DEGsMatrix_rna <- DEGsMatrix_rna[,order(colnames(DEGsMatrix_rna))]
    
    #batch effect treatment
    #matrix <- rbind(DEGsMatrix_mic, DEGsMatrix_rna)
    #labels <- c(labels_mic,labels_rna)
    #output$dim1<- renderText(dim(DEGsMatrix_rna))
    #output$dim2<- renderText(dim(DEGsMatrix_mic))
    #matrix <- batchEffectRemoval(as.matrix(t(matrix)), labels, method = "sva") #BATCH EFFECT REMOVAL
    #DEGsMatrix_mic <- matrix[,1:814]
    #DEGsMatrix_rna <- matrix[,815:dim(matrix)[2]]
    
    w$hide()
    
  #values$DEGsMatrix_rna_batch <- DEGsMatrix_rna
    
    #DEGs extraction
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running DEGs 10-folds CV..."),
                                        style="color:white;")))
    w$show
    set.seed(2)
    
    DEGs_CV <- DEGsExtraction(as.matrix(t(DEGsMatrix_mic)), as.factor(labels_mic), lfc=as.numeric(input$LFC), cov = as.numeric(input$COV), pvalue = as.numeric(input$pvalue), number = Inf, CV = TRUE, numFolds=10)
    DEGs <- DEGs_CV$Common_DEGs
    output$numberDEGs<- renderText(paste0( length(DEGs), " DEGs were extracted "))
    #huella_limma_mrmr_tri_mic_1a_2 <- featureSelection(t(expressionMatrix_tri_mic_1_2), labels_tri_mic_1, DEGs_2, mode = 'mrmr')
    w$hide()
    # mRMR method
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running mRMR algorithm..."),
                                        style="color:white;")))
    w$show()
    mrmrRanking <- featureSelection(as.matrix(DEGsMatrix_mic), as.factor(labels_mic), DEGs,
                                    mode = "mrmr")
    mrmrRanking <- names(mrmrRanking)
    w$hide()
    
  values$ranking <- mrmrRanking
  
    
  # Ranking tables
    
  output$genes_mrmr <- renderTable({
    return(mrmrRanking)
  }, colnames = FALSE)
  

  }) # Close button

  # Server of tab: Model training ------
  
  w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                 span(br(), br(), br(), h4(""),
                                 style="color:white;")))  
  
  observeEvent(input$boton_model_training, {
    
    w2$show()
    
    
    labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
    DEGsMatrix_rna <- as.matrix(read.csv(file = input$file_DEGsMatrix_rna$datapath))
    
    set.seed(333)
    indices_rna <- createDataPartition(labels_rna, p = .80, list = FALSE,  times = 1)
     
    particion.entrenamiento.rna <- DEGsMatrix_rna[indices_rna,]
    particion.test <- DEGsMatrix_rna[-indices_rna,]
    
    # Labels
    labels_train_rna <- labels_rna[indices_rna]
    labels_test_rna  <- labels_rna[-indices_rna]
    

    w2$hide()
    

    ranking <- values$ranking
  
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Training kNN algorithm..."),
                                           style="color:white;")))  
    w3$show()
    set.seed(333)
    results_cv <- knn_trn(as.matrix(particion.entrenamiento.rna), labels_train_rna, ranking,
                                numFold = as.numeric(input$number_folds))
    values$optimalkNN_train <- results_cv$bestK
      
    w3$hide()
    
    output$optimal_knn <- renderText(paste0("\nOptimal number of neighbours k = ", results_cv$bestK))
    
    output$train <- renderPlot({
      plot(results_cv$accuracyInfo$meanAccuracy[1:20], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.79,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
      lines(results_cv$sensitivityInfo$meanSensitivity[1:20], col='blue', lwd=2, lty=2)
      lines(results_cv$specificityInfo$meanSpecificity[1:20], col='#FF8B00', lwd=2, lty=4)
      lines(results_cv$F1Info$meanF1[1:20], col='red', lwd=2, lty=4)
      legend(x=15.9 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)
      
      
    })
    
  }) 
    
    w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4(""),
                                         style="color:white;")))  
    observeEvent(input$boton_boxplot, {
      w2$show()
      
      labels_mic <- as.vector(t(read.csv2(file = input$file_labels_mic$datapath)))
      DEGsMatrix_mic <- as.matrix(read.csv(file = input$file_DEGsMatrix_mic$datapath))
      labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
      DEGsMatrix_rna <- as.matrix(read.csv(file = input$file_DEGsMatrix_rna$datapath))
    
      w2$hide()
    
      ranking <- values$ranking
      
      w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4(""),
                                           style="color:white;")))
      
      w2$show()
      
    
      genes_tri <- ranking[1:as.numeric(input$nboxplot)]
      labels_tri <- labels_mic
      for (i in 1:length(labels_tri)){
      if (labels_tri[i]=='Control'){
        labels_tri[i] <- 'Control_mic'
      } else if (labels_tri[i]=='SCC') {
        labels_tri[i] <- 'SCC_mic'
      } else {
        labels_tri[i] <- 'ACC_mic'
      } 
    }
    
      labels_tri_1 <- labels_rna
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
      box_mic_tri <- DEGsMatrix_mic[,which(colnames(DEGsMatrix_mic) %in% genes_tri)]
      box_rna_tri <- DEGsMatrix_rna[,which(colnames(DEGsMatrix_rna) %in% genes_tri)]
      box_matrix_tri <- rbind(box_mic_tri, box_rna_tri)
      labels_box_tri_a <- labels_box_tri[-c(815,849,896)]
      labels_box_tri_a <- c(labels_box_tri_a[1:814],labels_box_tri[c(896,815,849)], labels_box_tri_a[815:length(labels_box_tri_a)])
      box_matrix_tri_a <- box_matrix_tri[-c(896,815,849),]
      box_matrix_tri_a <- rbind(box_matrix_tri_a[1:814,],box_matrix_tri[c(896,815,849),], box_matrix_tri_a[815:dim(box_matrix_tri_a)[1],])
    
    w2$hide()
    
    output$boxplot <- renderPlot({
      dataPlot(t(as.matrix(box_matrix_tri_a)),labels_box_tri_a,mode = "genesBoxplot",toPNG = FALSE, colours = c("tomato", "chartreuse", "cyan", "red",'green4','blue'))
    })
    })
    
    
  
  
  
  # Server of tab: Model validation  ------
  
  w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Validating model..."),
                                  style="color:white;")))  
  
  observeEvent(input$boton_model_validation, {
    
    w3$show()
    
    labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
    DEGsMatrix_rna <- as.matrix(read.csv(file = input$file_DEGsMatrix_rna$datapath))
    
    set.seed(333)
    indices_rna <- createDataPartition(labels_rna, p = input$porcentaje_entrenamiento / 100, list = FALSE,  times = 1)
    particion.entrenamiento.rna <- DEGsMatrix_rna[indices_rna,]
    particion.test <- DEGsMatrix_rna[-indices_rna,]
    
    # Labels
    labels_train_rna <- labels_rna[indices_rna]
    labels_test_rna  <- labels_rna[-indices_rna]
    w3$hide()
    
    
    ranking <- values$ranking
    
    
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Validating kNN algorithm..."),
                                           style="color:white;")))  
    w3$show()
    set.seed(333)
    results_validation <- knn_test(train = particion.entrenamiento.rna, labels_train_rna,
                                     test = particion.test, labels_test_rna,
                                     ranking, bestK = values$optimalkNN_train)
    w3$hide()
    
    
    output$results_validation <- renderPlot({
      tabla <- results_validation$cfMats[[input$numero_genes_validation]]$table
      plotConfMatrix(tabla)
    })
    
  })
  
  w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Calculating AUC..."),
                                       style="color:white;")))  
  
  
  #AUC
  observeEvent(input$boton_calculate_AUC, {
    
    w4$show()
    
    labels_rna <- as.vector(t(read.csv2(file = input$file_labels_rna$datapath)))
    DEGsMatrix_rna <- as.matrix(read.csv(file = input$file_DEGsMatrix_rna$datapath))
    
    set.seed(333)
    indices_rna <- createDataPartition(labels_rna, p = input$porcentaje_entrenamiento / 100, list = FALSE,  times = 1)
    particion.entrenamiento.rna <- DEGsMatrix_rna[indices_rna,]
    particion.test <- DEGsMatrix_rna[-indices_rna,]
    
    # Labels
    labels_train_rna <- labels_rna[indices_rna]
    labels_test_rna  <- labels_rna[-indices_rna]
    w4$hide()
    
    
    ranking <- values$ranking
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Calculating AUC..."),
                                         style="color:white;")))
    
    w4$show()
    set.seed(333)
    response <- as.factor(labels_test_rna)
    aucs <- rep(NA, length(levels(response))) # store AUCs
    legendLabels <- as.character()
    colours <- c('red','blue','green')
    a <- list()
    
    for (i in seq_along(levels(response))) {
      cur.class <- levels(response)[i]
      binaryTraining.labels <- as.factor(labels_train_rna == cur.class)
      binaryTest.labels <- as.factor(labels_test_rna == cur.class)
      
      tri_knn_cv_mrmr_results <- knn_trn(as.matrix(particion.entrenamiento.rna), binaryTraining.labels, ranking)
      
      tri_knn_test_mrmr_results <- knn_test(as.matrix(particion.entrenamiento.rna), binaryTraining.labels, as.matrix(particion.test), binaryTest.labels, ranking, bestK = tri_knn_cv_mrmr_results$bestK)
      
      score <- tri_knn_test_mrmr_results$predictions[[as.numeric(input$numero_genes_validation)]]
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
      a[[i]] <- roc.x
      a[[i+3]] <- roc.y
      # store AUC
      auc <- performance(pred, "auc")
      auc <- unlist(slot(auc, "y.values"))
      aucs[i] <- auc
      legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
    }
    
    output$meanAUC <- renderText((paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))) 
    
    w4$hide()
    
    output$AUC <- renderPlot({
      plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Sensitivity", xlab="1 - Specificity", bty='n', cex.lab=1.3, cex.axis=1.3)
      lines(x=c(0,1), c(0,1))
      lines(a[[4]] ~ a[[1]], col = colours[1], lwd = 2)
      lines(a[[5]] ~ a[[2]], col = colours[2], lwd = 2)
      lines(a[[6]] ~ a[[3]], col = colours[3], lwd = 2)
      legend(x=0.48 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)
      
      
    })
    
    
  })
  
}

shinyApp(ui, server)