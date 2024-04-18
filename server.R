#Installing packages
library(data.table)
library(shiny)
library(shinyjs)
library(shinyBS)
library(httr)
library(DT)
library(RCurl)
library(BiocManager)
library(GEOquery)
library(limma)
library(oligo)
library(annotate)
library(ggplot2)

source('helpers.R')

#setwd('../GSE')

shinyServer(function(input, output) {
  shinyjs:: useShinyjs()  # Include shinyjs
  shinyjs::showElement("wait_msg5")
  
  #Deleting files from working directory
  print(paste0("working directory is ", getwd()))
  currentDir<- "./"
  ProjectDir<- dir (path = currentDir, pattern ="^Project", full.names = TRUE)
  #unlink(paste0(ProjectDir,"*"), recursive = TRUE)
  unlink(file.path(getwd(),c("*.CEL", "*.gz", "*.txt.gz","*.txt","*.jpg", "*.png", "*.html", "*.tar", "*.zip"))) 
  print("*.CEL, *.tar, *.gz, *.jpg,*.png, *.html, *.txt.gz, *.txt, and *.zip are gone")
  print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
  wwwDir<- "./www"    
  unlink(file.path(wwwDir,c ("*.html")))
  unlink(file.path(wwwDir,c ("*.zip")))
  
  #changing the file size limit
  options(shiny.maxRequestSize = 50*2048^2) #25MB per file is allowed
  
  #preventing the shiny app from being grayed out
  autoInvalidate <- reactiveTimer(10000)
  observe({
    autoInvalidate()
    cat(".")
  })
  
  #View user's guide
  output$users_guide_msg <- renderUI({tagList(tags$p(tags$b(tags$h2("User's guide"))))})
  output$users_guide_view <- renderUI({ tags$iframe(style="height:800px; width:100%", src="Users_guide.pdf")})
  
  #Links to GitHub and DockerHub
  GitHub <- a ("GitHub link", href="https://github.com/annakatsiki/dexplore", target='_blank')
  output$GitHub <- renderUI({tagList(tags$p(GitHub))})
  DockerHub <- a ("DockerHub link", href="https://hub.docker.com/r/akatsiki/dexplore", target='_blank')
  output$DockerHub <- renderUI({tagList(tags$p(DockerHub))})
  
  #About us
  output$about <- renderUI({tagList(tags$p(tags$b("DExplore: an online tool for detecting differentially expressed genes from mRNA microarray experiments")),
                                    tags$p(tags$u("Anna D. Katsiki"), tags$sup("1"),", Pantelis E. Karatzas", tags$sup("2"), ", Alexandros A. Georgakilas", tags$sup("3"), "and Constantinos E. Vorgias", tags$sup("1")),
                                    tags$p(tags$sup("1"), "National and Kapodistrian University of Athens, Faculty of Biology, Department of Biochemistry and Molecular Biology, Panepistimiopolis Athens, Greece"),
                                    tags$p(tags$sup("2"),"National Technical University of Athens (NTUA), School of Chemical Engineering, Department of Process Analysis and Plant Design, Unit of Process Control & Informatics, Zografou Campus, Athens, Greece"),
                                    tags$p(tags$sup("3"),"National Technical University of Athens (NTUA), School of Applied Mathematical and Physical Sciences, Physics Department, DNA Damage Laboratory, Zografou Campus, Athens, Greece"))})
  
  #check if all mandatory fields have a value
  observe({
    mandatoryFilled <- vapply(fieldsMandatory, function(x) {!is.null(input[[x]]) && input[[x]] != ""}, logical(1))
    mandatoryFilled <- all(mandatoryFilled)
    mandatory_upload_Filled <- vapply(fieldsMandatory_up, function(x) {!is.null(input[[x]])},logical(1))
    mandatory_upload_Filled <- all(mandatory_upload_Filled)
    
    #enable/disable the submit button
    shinyjs::toggleState(id = "go", condition = mandatoryFilled)
    shinyjs::toggleState(id = "up_go", condition = mandatory_upload_Filled)
  })
  
  #after submission of a GSE accession number
  observeEvent (input$go, {
    print("submit")
    shinyjs::hideElement ("firstInput1")
    shinyjs::showElement ("wait_msg4")
    shinyjs::hideElement ("firstInput2")
    shinyjs::showElement ("firstRun")
    shinyjs::disable ("GSE")
    shinyjs::disable ("go")
    shinyjs::disable("up_go")
    shinyjs::showElement ("submit_msg")
    shinyjs::hideElement ("error")
    shinyjs::showElement ("wait_msg")
    shinyjs::showElement("wait_msg5")
    
    #GSE accession number validation and downloading
    accession<- input$GSE
    print(paste0("The length of input$GSE is ", nchar(accession)))
    ftp<- createFtp(accession)
    
    if (validation(ftp)) {
      link <-createLink(input$GSE)
      gse<- createGSE(input$GSE)
      output$url <- renderUI({tagList(tags$p("URL link:", link))})
      shinyjs::showElement ("goToNextTab")
      
      #Choosing platform
      plat_choices <- c()
      choose_platform <- FALSE
      shinyjs::showElement("choosing_platform")
      shinyjs::hideElement ("wait_msg4")
      if (length(gse)>0) {
        choose_platform <- TRUE
        for (i in 1:length(gse)) { plat_choices <- c(plat_choices, annotation(gse[[i]])) } #class(plat_choices) >>> character
      }
      plat_df <- as.data.frame(plat_choices)
      for (i in 1:length(plat_choices)) { row.names(plat_df)[i]<- plat_choices[i] }
      output$choosePlatform<- renderUI({
        selectInput ( inputId = "choosePlatform",
                      label = "Choose the platform for the analysis",
                      choices = row.names(plat_df),
                      selected = row.names(plat_df)[1])
      })
      observeEvent(input$submit_platform, {
        shinyjs::disable("submit_platform")
        cross<- FALSE
        if (length(gse)>1) {
          cross <- TRUE
          idx <- match(input$choosePlatform, plat_choices)
          cat("Selected platform: ")
          cat(idx)
          cat("\n")
        }
        if (length(gse)==1) {
          idx<-1
        }
        
        phenodata<- pData(gse[[idx]])            
        df<- data.frame(sample= phenodata$geo_accession, treatment=phenodata$title)
        gse_names <- paste0("GSE",input$GSE)
        dir <- paste0("./",gse_names)
        
        shinyjs::showElement ("data_table")
        
        #pdata, pDataRel, newpdata
        for (i in 1:ncol(phenodata)) {
          name <- paste(colnames(phenodata[i]))
          if (name == "source_name_ch1") {
            colnames(phenodata)[i] <- "source"
            name <- "source"
          }
          col <- phenodata[i]
          grp <- grep("characteristics", name)
          if (i == 1) {
            if (length(grp) == 0) {
              col <- sub("^", paste(name, ": ", sep = ""), col[, c(name)])
            }
            pdata <- as.data.frame(col)
            colnames(pdata)[i] <- name
          }
          else {
            if (length(grp) == 0) {
              col <- sub("^", paste(name, ": ", sep = ""), col[, c(name)])
            }
            pdata <- cbind(pdata, col)
            colnames(pdata)[i] <- name
          }
        }
        pdata <- as.data.frame(pdata)
        pdataRel <- pdata #relevant P data to be formed
        nameVec = c()
        cat("************************************************* \n")
        show(head(pdata))
        for (i in names(pdata)) {
          col <- pdata[, c(i)]
          title <- strsplit(paste(col[1]), ":")[[1]][1]
          numGroups <- length(unique(col))
          if ((numGroups < length(col) && numGroups > 1 && !(title %in% nameVec)) || title == "title") {
            nameVec <- c(nameVec, title)
            cat(
              "Column name: ",
              title,
              "\n",
              "Number of different groups: ",
              numGroups,
              "\n",
              sep = ""
            )
            if (numGroups < 21) {
              cat("Here is a list of groups and their frequencies: ",
                  "\n",
                  "\n")
            }
            if (numGroups >= 21) {
              cat("Too many groups to display. Here is a sample of the groups: ",
                  "\n")
            }
            freq_table <- as.data.frame(table(col))
            freq_table <- freq_table[order(-freq_table$Freq), c(1, 2)]
            colnames(freq_table) <- c(title, "Freq")
            for (j in 1:21) {
              if (j > numGroups) {
                break
              }
              group <- paste(freq_table[j, 1])
              cat(group, ": ", paste(freq_table[j, 2]), "\n")
            }
            cat("===============================================",
                "\n")
          }
          else {
            pdataRel[, c(i)] <- NULL
          }
        }
        colnames(pdataRel) <- nameVec
        pdata <- as.data.frame(pdataRel, stringsAsFactors = FALSE)
        newpdata <- sapply(pdata, as.character)
        show(head(pdata))
        if (ncol(pdata) == 0) {
          pdata =  as.data.frame(matrix(
            data = "",
            nrow = nrow(pdata),
            ncol = 1
          ))
          newpdata = pdata
          colnames(pdata[1]) = "No characteristics found"
          colnames(newpdata[1]) = "No characteristics found"
        } 
        else {
          for (i in 1:(ncol(pdata))) {
            for (j in 1:(nrow(pdata))) {cont <- as.character(strsplit(paste(pdata[j, i]), ":")[[1]]) 
            newpdata[j, i] <- cont[length(cont)]}
          }
        }
        rownames(newpdata) <- phenodata$geo_accession
        newpdata <- gsub("^ ", "", newpdata)
        output$df<- renderDataTable(newpdata, server = TRUE)

        #Choosing the comparison to be performed
        shinyjs::showElement ("comparison")
        output$compBut <- renderUI ({ 
          radioButtons ( inputId = "comp2",
                         label = "Select the criterion of the comparison",
                         choices = (colnames(newpdata)))    })
        shinyjs::enable("submit_comp")
        
        observeEvent(input$submit_comp, {
          shinyjs::disable("submit_comp")
          shinyjs::showElement ("ctr")
          output$comp_ctr <- renderUI ({ 
            radioButtons ( inputId = "ctr",
                           label = "Select which samples should be treated as control",
                           choices = (unique(newpdata[,input$comp2]))
            )
          })
          shinyjs::enable("submit_ctr")
          
          observeEvent(input$submit_ctr, {
            shinyjs::disable("submit_ctr")
            shinyjs::showElement ("tr")
            output$comp_tr <- renderUI ({
              radioButtons ( inputId = "tr",
                             label = "Select the samples to which control are to be compared",
                             choices = (unique(newpdata[which(!(newpdata[,input$comp2] %in% input$ctr)),input$comp2]))
              ) 
            })
            shinyjs::enable("submit_tr")
            
            observeEvent(input$submit_tr, {
              shinyjs::disable("submit_tr")
              print("comp_submitted")
              index <- which(colnames(pdataRel)==input$comp2)
              print(paste0("index is ",index))
              shinyjs::disable("submit_comp")
              shinyjs::showElement ("btns")
              shinyjs::showElement("run")
              if ("replicate" %in% colnames(newpdata)) { }
              
              else {
                replicate<- matrix(data = integer(), nrow = dim(newpdata)[1], ncol = length(unique(newpdata[, index])))
                row.names(replicate) <- row.names(newpdata)
                colnames(replicate) <- unique(newpdata[, index])
                replicate[1,1]<-1
                for (i in 2:dim(newpdata)[1]){
                  c<- which(colnames(replicate)==newpdata[i, index], arr.ind = T)
                  if (! is.na(replicate[i-1,c])) {replicate[i,c]<- replicate[i-1,c]+1}
                  else {replicate[i,c]<-1}
                }
                replicate
                replicate<- as.integer(replicate)
                replicate<- replicate[!is.na(replicate)]
                newpdata <- cbind(newpdata, replicate= replicate)
              }
              
              ###Creating targets.txt
              targets<-data.frame(newpdata)
              
              #run the analysis
              observeEvent(input$anls,{
                print("running the analysis")
                shinyjs::disable ("changes_done")
                shinyjs::disable ("anls")
                shinyjs::showElement ("wait_msg2")
                shinyjs::hideElement ("firstRun")
                shinyjs::showElement ("wait_msg3")

                ### downloading and untarring RAW files
                getGEOSuppFiles(gse_names, makeDirectory = FALSE)
                tar_file <- paste0(gse_names,"_RAW.tar")
                untar(tar_file, exdir = getwd())
                cF<- data.frame(list.files(pattern = "CEL[.]gz$"), stringsAsFactors = FALSE)
                print(cF)
                i<-1
                while (i<= dim(cF)[1]) {
                  gunzip(cF[i,], remove = FALSE, overwrite = T)
                  print(i)
                  i = i+1
                }
                
                ###Reading targets.txt
                print(targets)
                
                ###Reading Affy
                lc<- list.celfiles(path = getwd(), full.names = TRUE)
                celFiles<- c()
                for (i in 1:length(phenodata$geo_accession)) {
                  celFiles<- c(celFiles, grep(phenodata$geo_accession[i], lc, value = T))
                  }
                Data<- read.celfiles(celFiles)
                print("Data is ready")
                Data_rma<-rma(Data)
                print("Data_rma is ready")
                e <- exprs(Data_rma)
                print("e is ready")
                print(dim(e))
                
                ###NonSpecific Filtering-pOverA

                num_ctr <- length(which(newpdata[,index] %in% input$ctr)) ###number of control samples
                num_tr <- length(which(newpdata[,index] %in% input$tr)) ###number of treated samples
                prop<- (num_tr/(num_tr+num_ctr))
                print(paste("num_ctr is", num_ctr, sep = " "))
                print(paste("num_tr is", num_tr, sep = " "))
                print(paste("prop is", prop, sep = " "))
                f2_Data_rma<-filter1(Data_rma,prop,log2(100))  ###prop=number of treated/sum of samples
                dim(f2_Data_rma)
                f2_e<- filter2(e, prop, log2(100))
                dim(f2_e)
                
                ###NonSpecific Filtering-shorth
                row.mean <- esApply(Data_rma,1,mean) 
                sh <- shorth(row.mean)
                #hist(e)
                #abline(v=sh, col="red")
                f3_Data_rma <- Data_rma[row.mean >=sh,]
                dim(f3_Data_rma)
  
                ###Filtering Results - Comparison
                v<- c(dim(Data_rma)[1], dim(f2_Data_rma)[1], dim(f3_Data_rma)[1])
                filt.matrix<- data.frame(features=v, filt.method = c("before filt", "pOverA", "shorth"))  
                print(filt.matrix)
                
                ###Linear Modelling
                cond<-pdataRel[,index]
                batch<-as.factor(targets$replicate)
                design<-model.matrix(~0 + cond + batch)
                print("the design matrix before:")
                print (design)
                x_index <- grep(input$ctr, colnames(design), ignore.case = T) ###
                y_index <- grep(input$tr, colnames(design), ignore.case = T) ###
                colnames(design)<-gsub("cond","",colnames(design))
                colnames(design)<-gsub("atch","",colnames(design))
                colnames(design)<- make.names(colnames(design))
                print("the design matrix after:")
                print (design)
                fit<- lmFit(f2_e, design)
                print(head(fit$coefficients))
                x<- colnames(design)[x_index] ###
                y<- colnames(design)[y_index] ###
                s<- paste0(x,y,"=",y,"-",x)
                contr<- do.call(makeContrasts, c(s,list(levels = design)))
                fit2<-contrasts.fit(fit, contr)
                fit3<-eBayes(fit2,trend=T, robust = T)
                top<-topTable(fit3, number = nrow(fit3), adjust.method = input$adjMet, sort="none")
                head(top)
                dim(top)

                ###Significant DE genes
                print(input$adjMet)
                print(input$adjPval)
                print(input$logFC)
                res<-decideTests(fit3,method = "separate", adjust.method = input$adjMet, p.value = input$adjPval, lfc = input$logFC)
                
                # PLOTS
                
                # Histogram

                  png("histogram.png")
                  par(font.lab=2, font.axis=2)
                  hist(top$adj.P.Val, col = "darkgrey", border = "white", 
                       breaks = 100, 
                       xlab = "P-adj", ylab = "Number of probes", 
                       main = paste0("GSE", input$GSE, " adjusted P-value distribution"),
                       xaxt= "n"
                  )
                  axis(side = 1, at = seq(0, 1, length.out = 5), 
                       labels = seq(0, 1, length.out = 5), tck= 0)
                  dev.off()

                  # Boxplot
                  
                  png("boxplot.png")
                  ctr_index<-which(newpdata[,index] %in% input$ctr)
                  tr_index<-which(newpdata[,index] %in% input$tr)
                  groups<- c("Control", "Treated")
                  test_f2_e<-f2_e
                  sample_tags <- rep("Unknown", ncol(test_f2_e))
                  for (i in ctr_index) {
                    sample_tags[i] <- "Control"
                    }
                  for (i in tr_index) {
                    sample_tags[i] <- "Treated"
                    }
                  others <- setdiff(colnames(test_f2_e), c("Control", "Treated"))
                  colnames(test_f2_e)<-sample_tags
                  unknown_v<-(which(colnames(test_f2_e)=="Unknown"))
                  cols<-ifelse(colnames(test_f2_e)=="Control", "lightgreen", 
                               ifelse(colnames(test_f2_e)=="Treated", "red", "lightgrey"))
                  xlabels<-rownames(targets)
                  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, bty="l", font.lab=2, font.axis=2)
                  boxplot(f2_e, ylab="Expression Value", xaxt="n", 
                                col=cols,
                                main = paste0("Boxplot of GSE", input$GSE, axes=0)
                                )
                  text(x = 1:length(xlabels), y = par("usr")[3] - 0.5, 
                       labels = xlabels, srt=45, font=2,
                       adj = c(1, 1), xpd = TRUE, cex = 1)
                  if (length(unknown_v)>0){
                    leg <- legend("topright", pch = 19, col = c("lightgreen", "red", "lightgray"),
                                  legend = c(input$ctr, input$tr, "other"),
                                  fill = c("lightgreen", "red", "lightgray"),
                                  plot = FALSE)
                    
                    legend(x = (leg$rect$left + leg$rect$w) * 0.95, y = leg$rect$top,
                           pch = 19, col = c("lightgreen", "red", "lightgray"),
                           legend = c(input$ctr, input$tr, "other"),
                           bg= "transparent", border = NA, 
                           horiz = FALSE, box.lwd = 0, box.col = "white")
                         } else{
                           leg <- legend("topright", pch = 19, col = c("lightgreen", "red"),
                                         legend = c(input$ctr, input$tr),
                                         fill = c("lightgreen", "red"),
                                         plot = FALSE)
                           
                           legend(x = (leg$rect$left + leg$rect$w) * 0.95, y = leg$rect$top,
                                  pch = 19, col = c("lightgreen", "red"),
                                  legend = c(input$ctr, input$tr),
                                  bg= "transparent", border = NA, 
                                  horiz = FALSE, box.lwd = 0, box.col = "white")
                         }
                  dev.off()

                ###Keeping DE genes to a data.frame
                  sign<-subset(top, abs(top$logFC)> input$logFC & top$adj.P.Val < input$adjPval)
                  dim(sign)
                  print(paste0("Before annotation: nrow(sign) is ", nrow(sign))) 
                  print("head of sign is: ")
                  print(head(sign))
                
                ###Message for zero DE genes
                if (nrow(sign)==0){
                  print("nrow of sign is 0 - no DEGs")
                  shinyjs::hideElement ("wait_msg3")
                  shinyjs::hideElement ("wait_msg2")
                  shinyjs::showElement ("msg")
                  shinyjs::showElement ("noDEGs")
                }
                
                else {
                  ###Annotation using annotate
                  dict <- read.table("./plat_dict_full.tsv",
                                     header = T, sep = "\t", colClasses = "character")
                  gse_acc <- paste0("GSE", input$GSE)
                  gse<- getGEO(gse_acc, destdir = getwd(), GSEMatrix = T, getGPL = F)
                  a <- annotation(gse[[idx]])
                  p<- grep(a, dict$Accession, fixed = T) #ignore.case = T, 
                  x<- dict$db[p]
                  y<-x[1]
                  set<- c("-", "?", "")
                  
                  if(!is.element(y,set)){ #checking if annotation file exists
                    
                    if(!require(y, character.only = T)) {
                      BiocManager::install(y, character.only = T, dependencies = T, ask=FALSE)
                      library(y, character.only = T)
                    }
                    
                    sign<- annotf(sign, y)
                    top<- annotf(top,y)
                    
                    # PCA
                    
                    tr_f2_e <- t(f2_e)
                    
                    # Create a vector indicating the treatment group for each sample
                    sample_groups <- newpdata[, index, drop = FALSE]
                    # Use the ctr_index and tr_index to subset the columns of tr_f2_e
                    ctr_index<-which(newpdata[,index] %in% input$ctr)
                    tr_index<-which(newpdata[,index] %in% input$tr)
                    ctr_samples <- tr_f2_e[ctr_index, ]
                    trt_samples <- tr_f2_e[tr_index, ]
                    # Combine the control and treatment samples
                    selected_samples <- rbind(ctr_samples, trt_samples)
                    # Subset the sample groups vector accordingly
                    selected_sample_groups <- c(sample_groups[ctr_index], sample_groups[tr_index])
                    # Perform PCA
                    pca_result <- prcomp(selected_samples)
                    #print(summary(pca_result))
                    
                    ### Scree plot
                    #Calculate variance explained as percentage and cumulative variance
                    variance_explained_percent <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
                    cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
                    
                    # Create data frame for scree plot
                    scree_data <- data.frame(
                      PC = 1:length(pca_result$sdev),
                      Variance_Explained = variance_explained_percent,
                      Cumulative_Variance = cumulative_variance
                    )
                    
                    # Plot scree plot
                    scree_plot <- ggplot(scree_data, aes(x = PC)) +
                      geom_line(aes(y = Variance_Explained, color = "Variance Explained")) +
                      geom_line(aes(y = Cumulative_Variance, color = "Cumulative Variance Explained")) +
                      scale_color_manual(values = c("black", "red")) +
                      labs(title = paste0("Scree Plot - GSE", input$GSE), x = "Principal Component", y = "Variance Explained (%)") +
                      theme_minimal() +  ### change for uploaded data
                      theme(
                        plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        panel.border = element_blank(),  # Remove border
                        axis.line = element_line(color = "black"),  # Add axis lines with black color
                        axis.text = element_text(color = "black"),  # Set axis text color to black
                        axis.title = element_text(color = "black"),  # Set axis title color to black
                        plot.title = element_text(color = "black"),  # Set plot title color to black
                        plot.margin = margin(1, 1, 1, 1, "cm")  # Set plot margins
                      ) +
                      theme(
                        plot.title = element_text(face = "bold", hjust = 0.5),  
                        axis.title = element_text(face = "bold"),  # Bold axis labels
                        axis.title.x = element_text(hjust = 0.5),  # Centered x-axis label
                        axis.title.y = element_text(hjust = 0.5)   # Centered y-axis label
                      )
                    # Save scree plot
                    ggsave("scree_plot.png", plot = scree_plot, width = 8, height = 6, units = "in", dpi = 300)
                    
                    ### PCA_scores_plot
                    # Extract the scores (PC1 and PC2) and combine with sample names 
                    pc_df <- data.frame(
                      PC1 = pca_result$x[,1],
                      PC2 = pca_result$x[,2],
                      PC3 = pca_result$x[,3], 
                      Sample_Group = ifelse(seq_along(rownames(pca_result$x)) %in% ctr_index, "Control", "Treated")
                      )
                    # Plot PCA scores with sample groups
                    pca_scores_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Sample_Group)) +
                      geom_point() +
                      labs(title = paste0("PCA Plot with Sample Groups - GSE", input$GSE), x = "PC1", y = "PC2") +  ### change for uploaded data
                      scale_color_manual(values = c("Control" = "blue", "Treated" = "red")) + ### group names
                      theme_minimal() +
                      theme(
                        plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        panel.border = element_blank(),  # Remove border
                        axis.line = element_line(color = "black"),  # Add axis lines with black color
                        axis.text = element_text(color = "black"),  # Set axis text color to black
                        axis.title = element_text(color = "black"),  # Set axis title color to black
                        plot.title = element_text(color = "black"),  # Set plot title color to black
                        plot.margin = margin(1, 1, 1, 1, "cm")  # Set plot margins
                      ) +
                      theme(
                        plot.title = element_text(face = "bold", hjust = 0.5),  
                        axis.title = element_text(face = "bold"),  # Bold axis labels
                        axis.title.x = element_text(hjust = 0.5),  # Centered x-axis label
                        axis.title.y = element_text(hjust = 0.5)   # Centered y-axis label
                        )
                    # Save PCA scores plot
                    ggsave("pca_scores_plot.png", plot = pca_scores_plot + 
                             theme(panel.background = element_rect(fill = "white")), 
                           width = 8, height = 6, units = "in", dpi = 300)

                    ### Biplot
                    # Extract PCA scores and loadings
                    pc_scores <- as.data.frame(pca_result$x)
                    pc_loadings <- as.data.frame(pca_result$rotation)
                    annotated_pc_loadings <- annotf(pc_loadings, y)
                    # Filter out NAs
                    pc_loadings_no_na <- no_annotf(annotated_pc_loadings)
                    # Filter out duplicates
                    pc_loadings_filtered <- unFilt(annotated_pc_loadings$gene_symbol, pc_loadings_no_na)
                    # Choose top genes contributing to PC1 and PC2
                    top_genes <- 10  # Number of top genes to select
                    top_genes_pc1 <- order(abs(pc_loadings_filtered$PC1), decreasing = TRUE)[1:top_genes]
                    top_genes_pc2 <- order(abs(pc_loadings_filtered$PC2), decreasing = TRUE)[1:top_genes]
                    
                    Sample_Group = ifelse(seq_along(rownames(pca_result$x)) %in% ctr_index, "Control", "Treated")
                    
                    # Extract gene symbols
                    top_genes_symbols_pc1 <- pc_loadings_filtered$gene_symbol[top_genes_pc1]
                    top_genes_symbols_pc2 <- pc_loadings_filtered$gene_symbol[top_genes_pc2]
                    # Calculate the percentage of variance explained for PC1 and PC2
                    variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)
                    # Plot PCA loadings with biplot
                    biplot <- ggplot() +
                      labs(title = paste0("Biplot: PCA Loadings - GSE", input$GSE), ### change for uploaded data
                           x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
                           y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +  
                           theme_minimal() +
                      theme(
                        plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        panel.border = element_blank(),  # Remove border
                        axis.line = element_line(color = "black"),  # Add axis lines with black color
                        axis.text = element_text(color = "black"),  # Set axis text color to black
                        axis.title = element_text(color = "black"),  # Set axis title color to black
                        plot.title = element_text(color = "black"),  # Set plot title color to black
                        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")  # Set plot margins
                      ) +
                      scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  
                      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
                      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
                      coord_fixed(ratio = 1)
                    # Add arrows for top genes contributing to PC1
                    biplot <- biplot +
                      geom_segment(
                        data = data.frame(
                          x = 0, y = 0,
                          xend = pc_loadings_filtered[top_genes_pc1, "PC1"],
                          yend = pc_loadings_filtered[top_genes_pc1, "PC2"],
                          gene_symbol = top_genes_symbols_pc1
                          ),
                        aes(x = 0, y = 0, xend = xend, yend = yend),
                        arrow = arrow(length = unit(0.3, "cm")),
                        color = "darkgreen"
                        ) +
                      geom_text(
                        data = data.frame(
                          x = pc_loadings_filtered[top_genes_pc1, "PC1"],
                          y = pc_loadings_filtered[top_genes_pc1, "PC2"],
                          gene_symbol = top_genes_symbols_pc1
                          ),
                        aes(x = x, y = y, label = gene_symbol),
                        vjust = 1.5,
                        hjust = -0.5,
                        color = "darkgreen",
                        size = 2
                        )
                    # Add arrows for top genes contributing to PC2
                    biplot <- biplot +
                      geom_segment(
                        data = data.frame(
                          x = 0, y = 0,
                          xend = pc_loadings_filtered[top_genes_pc2, "PC1"],
                          yend = pc_loadings_filtered[top_genes_pc2, "PC2"],
                          gene_symbol = top_genes_symbols_pc2
                          ),
                        aes(x = 0, y = 0, xend = xend, yend = yend),
                        arrow = arrow(length = unit(0.3, "cm")),
                        color = "darkorange"
                        ) +
                      geom_text(
                        data = data.frame(
                          x = pc_loadings_filtered[top_genes_pc2, "PC1"],
                          y = pc_loadings_filtered[top_genes_pc2, "PC2"],
                          gene_symbol = top_genes_symbols_pc2
                          ),
                        aes(x = x, y = y, label = gene_symbol),
                        vjust = 1.5,
                        hjust = -0.5,
                        color = "darkorange",
                        size  = 2
                        )
                    # Customize axis labels to be bold and centered
                    biplot <- biplot +
                      theme(axis.title.x = element_text(face = "bold", hjust = 0.5),
                      axis.title.y = element_text(face = "bold", hjust = 0.5),
                      plot.title = element_text(face = "bold", hjust = 0.5))  
                    # Save biplot
                    ggsave("biplot.png", plot = biplot + 
                             theme(panel.background = element_rect(fill = "white")), 
                           width = 15, height = 8, units = "in", dpi = 300)

                    #Volcano plot
                    
                    library(plotly)
                    
                    # Plot each dataset with its corresponding color
                    # Create separate datasets for each category
                    upregulated <- top[top$logFC > input$logFC & top$adj.P.Val < input$adjPval, ]
                    downregulated <- top[top$logFC < -input$logFC & top$adj.P.Val < input$adjPval, ]
                    not_significant <- top[!(top$logFC > input$logFC & top$adj.P.Val < input$adjPval) & !(top$logFC < -input$logFC & top$adj.P.Val < input$adjPval), ]
                    max_y <- max(-log10(c(upregulated$adj.P.Val, downregulated$adj.P.Val, not_significant$adj.P.Val)))
                    
                    # Create custom hover text including gene symbol, log2FC, and adjusted p-value
                    up_text <- paste("Gene Symbol: ", upregulated$gene_symbol, "<br>",
                                     "log2FC: ", format(round(upregulated$logFC, 2), nsmall = 2), "<br>",
                                     "Adjusted P-value: ", format(round(upregulated$adj.P.Val, 3), nsmall = 3))
                    
                    down_text <- paste("Gene Symbol: ", downregulated$gene_symbol, "<br>",
                                       "log2FC: ", format(round(downregulated$logFC, 2), nsmall = 2), "<br>",
                                       "Adjusted P-value: ", format(round(downregulated$adj.P.Val, 3), nsmall = 3))
                    
                    not_significant_text <- paste("Gene Symbol: ", not_significant$gene_symbol, "<br>",
                                                  "log2FC: ", format(round(not_significant$logFC, 2), nsmall = 2), "<br>",
                                                  "Adjusted P-value: ", format(round(not_significant$adj.P.Val, 3), nsmall = 3))
                    
                    # Create the plot with custom hover text

                   plot <- plot_ly() %>%
                      add_trace(data = upregulated, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                                marker = list(size = 4, color = 'red'), name = 'Upregulated', text = up_text, hoverinfo = "text") %>%
                      add_trace(data = downregulated, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                                marker = list(size = 4, color = 'blue'), name = 'Downregulated', text = down_text, hoverinfo = "text") %>%
                      add_trace(data = not_significant, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                                marker = list(size = 2, color = 'black'), name = 'Not significant', text = not_significant_text, hoverinfo = "text") %>%
                      layout(
                        title = list(text = paste0("Interactive Volcano Plot", input$GSE), y = 0.97),  
                        xaxis = list(title = "log2 Fold Change", showgrid = FALSE, zeroline = FALSE), 
                        yaxis = list(title = "-log10(Adjusted P-value)", showgrid = FALSE, linecolor= "black"), 
                        hovermode = "closest",
                        legend = list(
                          orientation = "v",
                          x = 0.9,
                          y = 0.2,
                          xanchor = "left",
                          yanchor = "top",
                          bgcolor = "rgba(0,0,0,0)",
                          bordercolor = "rgba(0,0,0,0)",  
                          borderwidth = 2
                        ),
                        margin = list(l = 50, r = 50, b = 50, t = 50),  
                        font = list(family = "Arial", size = 14, color = "black", weight = "bold"),
                        showlegend = TRUE,
                        yaxis2 = list(range = c(0, max_y + 1)), 
                        print("volcano is ready")
                      )
                    
                    htmlwidgets::saveWidget(plot, file = "volcanoplot.html")
                    
                    
                    # Heatmap
                   
                    library(RColorBrewer)
                    library(heatmaply)
                    
                    sel.top <- subset(top, top$adj.P.Val< input$adjPval & abs(top$logFC)> input$logFC)
                    #sel.top<- top[which(abs(top$logFC)> input$logFC & top$adj.P.Val < input$adjPval),]
                    
                    
                    sel.top<- annotf(sel.top, y)
                    sel.top<- unFilt(sel.top$gene_symbol, sel.top)
                    sel.top<- sel.top[is.na(sel.top$gene_symbol)==FALSE,]
                    sel.top<- sel.top[order(sel.top$logFC, decreasing = TRUE),]
                    f2_e.sel<- f2_e[rownames(sel.top), ]
                    
                    keys<-sel.top$probeID
                    values<-sel.top$gene_symbol
                    l<- list()
                    for(i in 1:length(keys)){
                      l[keys[i]]<- values[i]
                      }
                    rownames(f2_e.sel)<- unlist(l[rownames(f2_e.sel)], use.names = F)
                    yellowblackblue <- colorRampPalette(c("dodgerblue", "black", "gold"))(n = 100)
                    data<- matrix(data = f2_e.sel, nrow = nrow(f2_e.sel))
                    hover_text <- matrix(paste("AveExpr:", round(data, 4)), 
                                         nrow = nrow(data), ncol = ncol(data))
                    
                    #subsetting the expression value matrix to view only the ctr and treated samples
                    selcol<- c (which(newpdata[, index] %in% input$ctr),which(newpdata[, index] %in% input$tr))
                    heatdata <- f2_e.sel[,selcol] 
                    
                    heat<- heatmaply(heatdata, distfun = "pearson", hclust_method = NA, 
                              colors = yellowblackblue, row_dend_left = T,
                              margins = c(50,50,70,0),
                              main = paste0("GSE",input$GSE), scale = "row", Colv = F,
                              label_names = c("Gene", "Sample", "Value"),
                              custom_hovertext = hover_text,
                              file = "./heatmap.html")
                    #browseURL("./heatmap.html")

                    files_to_zip<- c("volcanoplot.html",
                                     "boxplot.png",
                                     "histogram.png",
                                     "scree_plot.png",
                                     "pca_scores_plot.png",
                                     "biplot.png",
                                     "heatmap.html"
                                     )
                    zip_filename<- paste0("GSE", input$GSE,"_plots.zip")
                    zip(zip_filename, files_to_zip)
                    
                    output$down_plots <- downloadHandler(
                      filename = function() {
                        zip_filename
                      },
                      content = function(file) {
                        path_to_zip <- zip_filename
                        file.copy(path_to_zip, file)
                      }
                    )
                    
                    print(paste0("nrow(sign) after annotf is ",nrow(sign)))
                    
                    ###Filtering for NA geneNames
                    sign<- sign[is.na(sign$gene_symbol)==FALSE,]
                    print(paste0("nrow(sign) after removing NAs is ",nrow(sign)))
                    
                    ###Filtering for unique geneNames
                    sign_un<- unFilt(sign$gene_symbol, sign)
                    
                    print(paste0("DEGs table after unFilt:",nrow(sign_un)))
                    
                    if (nrow(sign_un)==0){
                      ###Message for zero DE genes
                      print("nrow of sign_un is 0 - no DEGs")
                      shinyjs::hideElement ("wait_msg3")
                      shinyjs::hideElement ("wait_msg2")
                      shinyjs::showElement ("msg")
                      shinyjs::showElement ("noDEGs")
                    }
                    else {
                      ###Rendering and downloading DE genes list
                      print(paste0("If the number of rows of sign_un is >=1: ", nrow(sign_un)))
                      sign_un <- sign_un[order(sign_un$adj.P.Val),]
                      degs <- sign_un
                      degs_t<- sign_un
                      sign_un$gene_symbol <- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/gene?term=", sign_un$gene_symbol,
                                                                        "[Gene Name] AND ", as.character(unique(phenodata$organism_ch1)), 
                                                                        "[Organism]"),
                                                    "' target='_blank'>", sign_un$gene_symbol,"</a>")
                      shinyjs::showElement ("DEGs")
                      sign_un<- DT::datatable( data = sign_un,                       
                                               editable = FALSE, 
                                               colnames = c("probeID","gene symbol","logFC","AveExpr",
                                                            "t","P.Value","adj.P.Val","B"),
                                               rownames = FALSE,
                                               escape = FALSE,
                                               extensions = c('Scroller' , 'Buttons'),
                                               options = list (
                                                 autoWidth = FALSE,
                                                 dom = 'Blftrip',
                                                 deferRender = T,
                                                 scrollY = 750,
                                                 scroller = T,
                                                 buttons = list(list(extend = 'colvis', columns = c(2, 3, 4, 5, 6, 7)),
                                                                "copy", "print")
                                               )
                      )%>% formatSignif(c('logFC','AveExpr','t','P.Value','adj.P.Val','B'),3)
                      
                      gene_link <- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/gene?term=", degs$gene_symbol,
                                                              "[Gene Name] AND ", as.character(unique(phenodata$organism_ch1)), 
                                                              "[Organism]"),
                                          "' target='_blank'>", degs$gene_symbol,"</a>")
                      degs$gene_symbol <- gene_link
                      degs <- cbind(degs, url = gene_link)
                      
                      output$DEGs <- DT::renderDataTable(sign_un, server = FALSE)
                      output$down_csv <- downloadHandler(
                        filename = function (){ paste0("GSE",input$GSE,"_DE_genes.csv")},
                        content = function (temp) { write.csv(degs, temp, row.names = FALSE)})
                      output$down_tsv <- downloadHandler(
                        filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                        content = function (temp) { write.table(degs, temp, quote = FALSE, sep = "\t",row.names = FALSE) })

                      
                      shinyjs::hideElement ("wait_msg2")
                      shinyjs::showElement ("msg")
                      shinyjs::hideElement ("wait_msg3")
                      shinyjs::showElement ("download")
                      shinyjs::enable ("download")
                      shinyjs::enable ("changes_done")
                      shinyjs::enable ("submit_comp")
                      shinyjs::enable ("anls")
                      shinyjs::hideElement ("wait_msg5")
                      shinyjs::showElement ("WebGestaltR_but")

                      library(WebGestaltR)

                      
                      observeEvent(input$WebGestaltR,{
                        shinyjs::showElement("WebGestaltR_but")
                        shinyjs::disable("WebGestaltR")
                        shinyjs::showElement("choosing_organism")
                        output$chooseOrganism<- renderUI({
                          selectInput( inputId= "chooseOrganism",
                                       label = "Choose the organism for the analysis",
                                       choices = listOrganism()
                          )
                        })
                        shinyjs::enable("submit_Organism")
                        print("Choosing organism")
                        
                        observeEvent(input$submit_Organism,{
                          shinyjs::disable("submit_Organism")
                          shinyjs::showElement("choosing_RefSet")
                          organismWebGestaltR<- input$chooseOrganism
                          output$chooseRefSet<- renderUI({
                            selectInput( inputId = "chooseRefSet",
                                         label = "Choose the reference set for the analysis",
                                         choices = listReferenceSet(organism = organismWebGestaltR))})
                          print("Choosing RefSet")
                          
                          observeEvent(input$submit_RefSet,{
                            shinyjs::disable("submit_RefSet")
                            shinyjs::hideElement("wait_msg5")
                            shinyjs::showElement("wait_msg6")
                            shinyjs::hideElement("WebGestalt_down")
                            shinyjs::hideElement("WebGestalt_res")
                            shinyjs::hideElement("WebGestaltR_but")
                            shinyjs::hideElement("choosing_organism")
                            shinyjs::hideElement("choosing_RefSet")
                            
                            library(zip)
                            #}
                            print("Running WebGestaltR")
                            genesymbol <- degs_t$gene_symbol
                            genesymbol <- as.character(genesymbol)
                            projectName <- paste0("GSE", input$GSE,"_WebGestalt")
                            outputDir <- currentDir
                            WebGestaltR(enrichMethod = "ORA", organism = input$chooseOrganism,
                                        enrichDatabase = c("geneontology_Biological_Process", "geneontology_Molecular_Function",
                                                           "geneontology_Cellular_Component"),
                                        interestGene= genesymbol, interestGeneType="genesymbol", 
                                        collapseMethod="mean", referenceSet= input$chooseRefSet,
                                        minNum=10, maxNum=500, sigMethod="fdr", fdrMethod="BH", fdrThr=0.05, 
                                        reportNum=20, isOutput= T, outputDirectory=outputDir,
                                        projectName=projectName, hostName= "https://www.webgestalt.org" )
                            HTMLfile<- list.files(path = paste0(outputDir,"/Project_",projectName), pattern = ".html")
                            print(paste0("HTML file is ",HTMLfile))
                            projectPath<- paste0(outputDir,"/","Project_", projectName)
                            print(paste0("project path is ", projectPath))
                            file.copy(from = paste0(projectPath,"/",HTMLfile),
                                      to = paste0(outputDir,"/www"))
                            shinyjs::hideElement("wait_msg6")
                            shinyjs::showElement("WebGestalt_down")
                            shinyjs::showElement("WebGestalt_res")
                            output$WebGestalt <- renderUI({
                              tags$iframe(style="height:800px; width:100%", src=HTMLfile)
                            })
                            Webgestalt_files<- list.files(path = paste0(outputDir,"/Project_", projectName), 
                                                          recursive = TRUE, full.names = TRUE)
                            zip::zipr(zipfile = "WebGestalt.zip", files = Webgestalt_files, include_directories=TRUE)
                            
                            if (file.exists("WebGestalt.zip")) {print ("zip is ready!!!")}
                            # the zip file is in working directory, here C:/Users/Anna/Desktop/PhD/DExplore/GSE
                            else {print("Problem while zipping")}
                          })
                          output$down_WebGestalt <- downloadHandler(
                            filename = function() {paste0("Project_GSE", input$GSE, "_WebGestalt.zip")},
                            content = function(temp) {file.copy(from = paste0(getwd(),"/WebGestalt.zip"), to = temp)} )
                        }) #submit_Organism
                      }) #WebGestaltR
                    } #Rendering and downloading DE genes list
                  } # if (!is.element(y,set))
                  
                  else { ### sto if(!is.element(y,set))
                    shinyjs::showElement("no_annotation")
                    ###Rendering and downloading DE genes list WITHOUT ANNOTATION 
                    sign<- no_annotf(sign)
                    dim(sign)
                    sign2 <- sign[order(sign$adj.P.Val),]
                    degs2<- sign2
                    shinyjs::showElement ("DEGs")
                    sign2<- DT::datatable( data = sign2,                       
                                           editable = FALSE, 
                                           colnames = c("probeID","logFC","AveExpr","t","P.Value","adj.P.Val","B"),
                                           rownames = FALSE,
                                           extensions = c('Scroller' , 'Buttons'),
                                           options = list (
                                             autoWidth = FALSE,
                                             dom = 'Blftrip',
                                             deferRender = T,
                                             scrollY = 750,
                                             scroller = T,
                                             buttons = list(list(extend = 'colvis', columns = c( 1, 2, 3, 4, 5, 6)),
                                                            "copy", "print")
                                           )
                    )%>% formatSignif(c('logFC','AveExpr','t','P.Value','adj.P.Val','B'),3)
                    output$DEGs <- DT::renderDataTable(sign2, server = FALSE)
                    output$down_csv <- downloadHandler(
                      filename = function (){ paste0("GSE",input$GSE,"_DE_genes.csv")},
                      content = function (temp) { write.csv(degs2, temp, row.names = FALSE)})
                    output$down_tsv <- downloadHandler(
                      filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                      content = function (temp) { write.table(degs2, temp, quote = FALSE, sep = "\t",row.names = FALSE) })
                    output$down_pdf <- downloadHandler( 
                      filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                      content = function ( temp ) {
                        createPdf(degs2)
                        file.copy ("test.pdf", temp)},
                      contentType = "application/pdf")
                    
                    shinyjs::hideElement ("wait_msg2")
                    shinyjs::showElement ("msg")
                    shinyjs::hideElement ("wait_msg3")
                    shinyjs::showElement ("download")
                    shinyjs::enable ("download")
                    shinyjs::enable ("changes_done")
                    shinyjs::enable ("submit_comp")
                    shinyjs::enable ("anls")
                    shinyjs::hideElement("wait_msg5")
                    shinyjs::showElement("msg_no_annot")
                  }  #else { sto if(!is.element(y,set))
                  
                  #Deleting files from working directory
                  unlink(file.path(getwd(),c("*.CEL", "*.tar", "*.gz", "*.txt.gz","*.txt"))) 
                  print("*.CEL, *.tar, *.gz, *.txt.gz and *.txt are gone")
                  print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
                  
                } #annotation
              }) #anls
            })#tr
          }) #ctr
        }) #comp
      }) #platform
    } #validation
    
    else { #validation(ftp) is FALSE
      output$falseGSE <- renderText({paste("GSE",input$GSE, " is not a valid accesion number. Please, check again.", sep = "")}) 
      shinyjs::showElement ("falseGSE")
      print("validation==FALSE")
    }
    
    shinyjs::hideElement ("submit_msg")
    shinyjs::hideElement ("wait_msg")
    
  }) #go
  
  
  ###### ***** FOR MANUALLY UPLOADED DATA ***** ######     
  
  
  #using uploaded data
  observeEvent(input$up_go,{
    print("uploaded_submit")
    
    shinyjs::hideElement ("firstInput1")
    shinyjs::hideElement ("firstInput2")
    shinyjs::showElement ("firstRun")

    shinyjs::disable ("GSE")
    shinyjs::disable ("go")
    shinyjs::disable ("go_up")
    shinyjs::hideElement ("error")
    
    fileNames <- matrix(input$uploadFiles$name)
    empty_up <- create_uploaded_df(fileNames)
    df_up <- create_uploaded_DataTable(fileNames)  
    output$df_up<- DT::renderDataTable(df_up, server = TRUE)
    
    #shinyjs::showElement ("goToNextTab")
    shinyjs::showElement ("data_table_up")
    shinyjs::showElement ("data_sub_msg")
    shinyjs::showElement ("changes")
    
    #Editing Data Description form
    observeEvent(input$df_up_cell_edit,{
      info <- input$df_up_cell_edit
      i <- info$row
      j <- info$col+1
      v <- info$value
      empty_up[i,j]<<- DT::coerceValue(v,empty_up[i,j])
    }) #df_up_cell_edit
    
    #Saving changes - Creating targets.txt
    observeEvent(input$changes_done,{
      shinyjs::disable("changes_done")
      print("changes_done")
      
      ###comparisons
      if(length(unique(empty_up$treatment))>1)     { a<- empty_up$treatment }      else { a <-"" }
      if(length(unique(empty_up$duration))>1)      { b<- empty_up$duration }       else { b <-"" }
      if(length(unique(empty_up$concentration))>1) { c<- empty_up$concentration }  else { c <-"" }
      abc <- paste(a, b, c, sep = "-")
      if ( (length(unique(abc))>=2) & (length(unique(empty_up$replicate))>=2) ){
        df_comp <- create_df_comp(unique(abc))
        print(paste0("df_comp is ", class(df_comp)))
        print(df_comp)
        m_comp <- create_m_comp(df_comp)
        print(paste0("m_comp is ", class(m_comp)))
        print(m_comp)
        shinyjs::hideElement ("no_comparison")
        shinyjs::showElement ("comparison")
        output$compBut <- renderUI ({ 
          radioButtons ( inputId = "comp2",
                         label = "Select the comparison",
                         choices = (m_comp)
                         ) 
          
        }) #comparison
        shinyjs::enable("submit_comp")
      }
      else  {
        shinyjs::showElement ("no_comparison")
        output$comp <- renderText( c ("Please, fill in the form correctly.",
                                      "Columns 'treatment' and 'replicate' should be filled out!",
                                      "You may use information from GEO's webpage."))
      } #no comparison
      emptyDT<- DT::datatable(empty_up,
                              editable = TRUE,
                              rownames = FALSE,
                              selection = "none",
                              extensions = 'Buttons',
                              options = list(autoWidth = FALSE,
                                             dom='Bfrtip',
                                             buttons = list(list(extend = 'colvis', columns = c(2, 3)))
                              )) %>% formatStyle('sample',  
                                                 color='black', 
                                                 backgroundColor = '#D1C4E9', 
                                                 fontWeight = 'bold')
      output$df <- DT::renderDataTable(emptyDT)
      empty_up$sample <- fileNames[,1]
      targets<-data.frame(empty_up)

      observeEvent(input$submit_comp, {
        print("comp_submitted")
        print("input$comp2")
        print(input$comp2)
        print(class(input$comp2))
        index <- which(m_comp==input$comp2)
        print(paste("index:", class(index), index, sep=" "))
        shinyjs::disable("submit_comp")
        shinyjs::showElement ("btns")
        shinyjs::showElement("run")
        
        #run the analysis
        observeEvent(input$anls,{
          print("running the analysis")
          print(index)
          print(df_comp[index,])
          shinyjs::disable ("changes_done")
          shinyjs::disable ("anls")
          shinyjs::showElement ("wait_msg2")
          shinyjs::hideElement ("firstRun")
          shinyjs::showElement ("wait_msg3")

          ###Reading Affy
          Data<- read.celfiles(input$uploadFiles$datapath)
          Data_rma<-rma(Data)
          dim(Data_rma)
          e <- exprs(Data_rma)
          print(dim(e))
          
          ###NonSpecific Filtering-pOverA
          num_tr<-length(subset(abc, abc==df_comp[index,3])) ###number of treated samples
          num_ctrl<-length(subset(abc, abc==df_comp[index,1])) ###number of control samples
          prop <- (num_tr/(num_tr+num_ctrl))
          print(paste0("prop is ",prop))
          f2_Data_rma<-filter1(Data_rma,prop,log2(100))  ###prop= number of treated/sum of samples
          dim(f2_Data_rma)
          f2_e<- filter2(e, prop, log2(100))
          dim(f2_e)
          
          ###NonSpecific Filtering-shorth
          row.mean <- esApply(Data_rma,1,mean) 
          sh <- shorth(row.mean)
          #hist(e)
          #abline(v=sh, col="red")
          f3_Data_rma <- Data_rma[row.mean >=sh,]
          dim(f3_Data_rma)
          
          ###Filtering Results - Comparison
          v<- c(dim(Data_rma)[1], dim(f2_Data_rma)[1], dim(f3_Data_rma)[1])
          filt.matrix<- data.frame(features=v, filt.method = c("before filt", "pOverA", "shorth"))  
          print(filt.matrix)
          
          ###Linear Modelling
          abc<-gsub("-","_",abc)
          cond<-as.factor(abc)
          batch<-as.factor(targets$replicate)
          design<-model.matrix(~0 + cond + batch)
          print("the design matrix before:")
          print (design)
          colnames(design)<-gsub("cond","",colnames(design))
          colnames(design)<-gsub("atch","",colnames(design))
          print("the design matrix after:")
          print (design)
          fit<- lmFit(f2_e, design)
          print(head(fit$coefficients))
          x<- gsub("-","_",df_comp[index,1])
          y<- gsub("-","_",df_comp[index,3])
          s<- paste0(x,y,"=",y,"-",x)
          contr<- do.call(makeContrasts, c(s,list(levels = design)))
          fit2<-contrasts.fit(fit, contr)
          fit3<-eBayes(fit2,trend=T, robust = T)
          top<-topTable(fit3, number = nrow(fit3), adjust.method = input$adjMet, sort="none")
          head(top)
          dim(top)
          
          ###Significant DE genes
          print(input$adjMet)
          print(input$adjPval)
          print(input$logFC)
          res<-decideTests(fit3,method = "separate", adjust.method = input$adjMet, p.value = input$adjPval, lfc = input$logFC)
          
          # PLOTS

          # Histogram
          
          png("histogram.png")
          par(font.lab=2, font.axis=2)
          hist(top$adj.P.Val, col = "darkgrey", border = "white", 
               breaks = 100, 
               xlab = "P-adj", ylab = "Number of probes", 
               main = "Adjusted P-value distribution",
               xaxt= "n"
          )
          axis(side = 1, at = seq(0, 1, length.out = 5), 
               labels = seq(0, 1, length.out = 5), tck= 0)
          dev.off()
          
          # Boxplot
          cond1<- df_comp[index,1]
          cond2<-df_comp[index,3]
          if(length(unique(empty_up$treatment))>1)     { a<- empty_up$treatment }      else { a <-"" }
          if(length(unique(empty_up$duration))>1)      { b<- empty_up$duration }       else { b <-"" }
          if(length(unique(empty_up$concentration))>1) { c<- empty_up$concentration }  else { c <-"" }
          abc <- paste(a, b, c, sep = "-")
          
          png("boxplot.png", width = 800, height = 600)                  
          ctr_index<- which(abc==cond1, arr.ind = TRUE)
          tr_index<- which(abc==cond2, arr.ind = TRUE)
          
          groups<- c("Control", "Treated")
          test_f2_e<-f2_e
          sample_tags <- rep("Unknown", ncol(test_f2_e))
          for (i in ctr_index) {
            sample_tags[i] <- "Control"
          }
          for (i in tr_index) {
            sample_tags[i] <- "Treated"
          }
          
          others <- setdiff(colnames(test_f2_e), c("Control", "Treated"))
          colnames(test_f2_e)<-sample_tags
          unknown_v<-which(colnames(test_f2_e)=="Unknown")
          
          cols<-ifelse(colnames(test_f2_e)=="Control", "lightgreen", 
                       ifelse(colnames(test_f2_e)=="Treated", "red", "lightgrey"))
          xlabels<-targets$sample
          par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, bty="l", font.lab=2, font.axis=2)
          boxplot(f2_e, ylab="Expression Value", xaxt="n", 
                  col=cols,
                  main = "Boxplot" )
          text(x = 1:length(xlabels), y = par("usr")[3] - 0.5, 
               labels = xlabels, srt=45, font=2,
               adj = c(1, 1), xpd = TRUE, cex = 1)
          cond1_clean <- gsub("[^[:alnum:]]", "", cond1)
          cond2_clean <- gsub("[^[:alnum:]]", "", cond2)
          
          if (length(unknown_v)>0){
            leg <- legend("topright", pch = 19, col = c("lightgreen", "red", "lightgray"),
                          legend = c(cond1_clean, cond2_clean, "other"), 
                          fill = c("lightgreen", "red", "lightgray"),
                          plot = FALSE)
            
            legend(x = (leg$rect$left + leg$rect$w) * 0.95, y = leg$rect$top,
                   pch = 19, col = c("lightgreen", "red", "lightgray"),
                   legend = c(cond1_clean, cond2_clean, "other"), 
                   bg= "transparent", border = NA, 
                   horiz = FALSE, box.lwd = 0, box.col = "white")
          } else {
            leg <- legend("topright", pch = 19, col = c("lightgreen", "red"),
                          legend = c(cond1_clean, cond2_clean), 
                          fill = c("lightgreen", "red"),
                          plot = FALSE)
            
            legend(x = (leg$rect$left + leg$rect$w) * 0.95, y = leg$rect$top,
                   pch = 19, col = c("lightgreen", "red"),
                   legend = c(cond1_clean, cond2_clean), 
                   bg= "transparent", border = NA, 
                   horiz = FALSE, box.lwd = 0, box.col = "white")
          }					
          
          dev.off()
          
          ###Keeping DE genes to a data.frame
          sign<-  subset(top, top$adj.P.Val< input$adjPval & abs(top$logFC)> input$logFC) 
          
          #sign<-top[which(abs(top$logFC) > input$logFC & top$adj.P.Val < input$adjPval),]
          print(paste0("Before annotation: nrow(sign) is ", nrow(sign))) 
          print("head of sign is: ")
          print(head(sign))
          
          ###Message for zero DE genes
          if (nrow(sign)==0){
            print("nrow of sign is 0 - no DEGs")
            shinyjs::hideElement ("wait_msg3")
            shinyjs::hideElement ("wait_msg2")
            shinyjs::showElement ("msg")
            shinyjs::showElement ("noDEGs")
          } #if (nrow(sign)==0){}
          
          else { 
            ###Annotation using annotate
            dict<- read.table("./plat_dict_full.tsv", 
                              header = T, sep="\t", colClasses = "character")
            p<- grep(Data@annotation,dict$annotation, ignore.case = T, fixed = T)
            x<- dict$db[p]
            y<-x[1]
            org<- dict$Taxonomy[p]
            org<- org[1]
            set<- c("-", "?", "")
            if(!is.element(y,set)){ #checking if annotation file exists
              if(!require(y, character.only = T)) {
                BiocManager::install(y, character.only = T, dependencies = T, ask = FALSE)
                library(y, character.only = T)
              }
              sign<- annotf(sign, y)
              print(paste0("nrow(sign) after annotf is ",nrow(sign)))
              
              top<- annotf(top,y)
              
              # PCA
              cond1<- df_comp[index,1]
              cond2<-df_comp[index,3]
              if(length(unique(empty_up$treatment))>1)     { a<- empty_up$treatment }      else { a <-"" }
              if(length(unique(empty_up$duration))>1)      { b<- empty_up$duration }       else { b <-"" }
              if(length(unique(empty_up$concentration))>1) { c<- empty_up$concentration }  else { c <-"" }
              abc <- paste(a, b, c, sep = "-")
              
              ctr_index<- which(abc==cond1, arr.ind = TRUE)
              tr_index<- which(abc==cond2, arr.ind = TRUE)
              groups<- c("Control", "Treated")
              pca_f2_e<-f2_e
              sample_tags <- rep("Unknown", ncol(pca_f2_e))
              for (i in ctr_index) {
                sample_tags[i] <- "Control"
              }
              for (i in tr_index) {
                sample_tags[i] <- "Treated"
              }
              others <- setdiff(colnames(pca_f2_e), c("Control", "Treated"))
              colnames(pca_f2_e)<-sample_tags
              
              # Identify samples belonging to control and treated groups
              control_samples <- which(sample_tags == "Control")
              treated_samples <- which(sample_tags == "Treated")
              
              # Subset the data to include only control and treated samples
              selected_samples <- pca_f2_e[, c(control_samples, treated_samples)]
              
              # Perform PCA
              pca_result <- prcomp(t(selected_samples))
              print(summary(pca_result))
              
              ### Scree plot
              #Calculate variance explained as percentage and cumulative variance
              variance_explained_percent <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
              cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
              
              # Create data frame for scree plot
              scree_data <- data.frame(
                PC = 1:length(pca_result$sdev),
                Variance_Explained = variance_explained_percent,
                Cumulative_Variance = cumulative_variance
              )
              
              # Plot scree plot
              scree_plot <- ggplot(scree_data, aes(x = PC)) +
                geom_line(aes(y = Variance_Explained, color = "Variance Explained")) +
                geom_line(aes(y = Cumulative_Variance, color = "Cumulative Variance Explained")) +
                scale_color_manual(values = c("black", "red")) +
                labs(title = paste0("Scree Plot"), x = "Principal Component", y = "Variance Explained (%)") +
                theme_minimal() +  
                theme(
                  plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                  panel.grid.major = element_blank(),  # Remove major grid lines
                  panel.grid.minor = element_blank(),  # Remove minor grid lines
                  panel.border = element_blank(),  # Remove border
                  axis.line = element_line(color = "black"),  # Add axis lines with black color
                  axis.text = element_text(color = "black"),  # Set axis text color to black
                  axis.title = element_text(color = "black"),  # Set axis title color to black
                  plot.title = element_text(color = "black"),  # Set plot title color to black
                  plot.margin = margin(1, 1, 1, 1, "cm")  # Set plot margins
                ) +
                theme(
                  plot.title = element_text(face = "bold", hjust = 0.5),  
                  axis.title = element_text(face = "bold"),  # Bold axis labels
                  axis.title.x = element_text(hjust = 0.5),  # Centered x-axis label
                  axis.title.y = element_text(hjust = 0.5)   # Centered y-axis label
                )
              # Save scree plot
              ggsave("scree_plot.png", plot = scree_plot, width = 8, height = 6, units = "in", dpi = 300)
              
              ### PCA_scores_plot
              # Extract the scores (PC1 and PC2) and combine with sample names 
              pc_df <- data.frame(
                PC1 = pca_result$x[,1],
                PC2 = pca_result$x[,2],
                PC3 = pca_result$x[,3], 
                Sample_Group = ifelse(seq_along(rownames(pca_result$x)) %in% ctr_index, "Control", "Treated")
              )
              # Plot PCA scores with sample groups
              pca_scores_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Sample_Group)) +
                geom_point() +
                labs(title = paste0("PCA Plot with Sample Groups"), x = "PC1", y = "PC2") +  
                scale_color_manual(values = c("Control" = "blue", "Treated" = "red")) + ### group names
                theme_minimal() +
                theme(
                  plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                  panel.grid.major = element_blank(),  # Remove major grid lines
                  panel.grid.minor = element_blank(),  # Remove minor grid lines
                  panel.border = element_blank(),  # Remove border
                  axis.line = element_line(color = "black"),  # Add axis lines with black color
                  axis.text = element_text(color = "black"),  # Set axis text color to black
                  axis.title = element_text(color = "black"),  # Set axis title color to black
                  plot.title = element_text(color = "black"),  # Set plot title color to black
                  plot.margin = margin(1, 1, 1, 1, "cm")  # Set plot margins
                ) +
                theme(
                  plot.title = element_text(face = "bold", hjust = 0.5),  
                  axis.title = element_text(face = "bold"),  # Bold axis labels
                  axis.title.x = element_text(hjust = 0.5),  # Centered x-axis label
                  axis.title.y = element_text(hjust = 0.5)   # Centered y-axis label
                )
              # Save PCA scores plot
              ggsave("pca_scores_plot.png", plot = pca_scores_plot + 
                       theme(panel.background = element_rect(fill = "white")), 
                     width = 8, height = 6, units = "in", dpi = 300)
              
              ### Biplot
              # Extract PCA scores and loadings
              pc_scores <- as.data.frame(pca_result$x)
              pc_loadings <- as.data.frame(pca_result$rotation)
              annotated_pc_loadings <- annotf(pc_loadings, y)
              # Filter out NAs
              pc_loadings_no_na <- no_annotf(annotated_pc_loadings)
              # Filter out duplicates
              pc_loadings_filtered <- unFilt(annotated_pc_loadings$gene_symbol, pc_loadings_no_na)
              # Choose top genes contributing to PC1 and PC2
              top_genes <- 10  # Number of top genes to select
              top_genes_pc1 <- order(abs(pc_loadings_filtered$PC1), decreasing = TRUE)[1:top_genes]
              top_genes_pc2 <- order(abs(pc_loadings_filtered$PC2), decreasing = TRUE)[1:top_genes]
              
              Sample_Group = ifelse(seq_along(rownames(pca_result$x)) %in% ctr_index, "Control", "Treated")
              
              # Extract gene symbols
              top_genes_symbols_pc1 <- pc_loadings_filtered$gene_symbol[top_genes_pc1]
              top_genes_symbols_pc2 <- pc_loadings_filtered$gene_symbol[top_genes_pc2]
              # Calculate the percentage of variance explained for PC1 and PC2
              variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)
              # Plot PCA loadings with biplot
              biplot <- ggplot() +
                labs(title = paste0("Biplot: PCA Loadings"), 
                     x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
                     y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +  
                theme_minimal() +
                theme(
                  plot.background = element_rect(fill = "white", color = NA),  # Set entire plot background to white
                  panel.grid.major = element_blank(),  # Remove major grid lines
                  panel.grid.minor = element_blank(),  # Remove minor grid lines
                  panel.border = element_blank(),  # Remove border
                  axis.line = element_line(color = "black"),  # Add axis lines with black color
                  axis.text = element_text(color = "black"),  # Set axis text color to black
                  axis.title = element_text(color = "black"),  # Set axis title color to black
                  plot.title = element_text(color = "black"),  # Set plot title color to black
                  plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")  # Set plot margins
                ) +
                scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  
                scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
                theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                      plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
                coord_fixed(ratio = 1)
              # Add arrows for top genes contributing to PC1
              biplot <- biplot +
                geom_segment(
                  data = data.frame(
                    x = 0, y = 0,
                    xend = pc_loadings_filtered[top_genes_pc1, "PC1"],
                    yend = pc_loadings_filtered[top_genes_pc1, "PC2"],
                    gene_symbol = top_genes_symbols_pc1
                  ),
                  aes(x = 0, y = 0, xend = xend, yend = yend),
                  arrow = arrow(length = unit(0.3, "cm")),
                  color = "darkgreen"
                ) +
                geom_text(
                  data = data.frame(
                    x = pc_loadings_filtered[top_genes_pc1, "PC1"],
                    y = pc_loadings_filtered[top_genes_pc1, "PC2"],
                    gene_symbol = top_genes_symbols_pc1
                  ),
                  aes(x = x, y = y, label = gene_symbol),
                  vjust = 1.5,
                  hjust = -0.5,
                  color = "darkgreen",
                  size = 2
                )
              # Add arrows for top genes contributing to PC2
              biplot <- biplot +
                geom_segment(
                  data = data.frame(
                    x = 0, y = 0,
                    xend = pc_loadings_filtered[top_genes_pc2, "PC1"],
                    yend = pc_loadings_filtered[top_genes_pc2, "PC2"],
                    gene_symbol = top_genes_symbols_pc2
                  ),
                  aes(x = 0, y = 0, xend = xend, yend = yend),
                  arrow = arrow(length = unit(0.3, "cm")),
                  color = "darkorange"
                ) +
                geom_text(
                  data = data.frame(
                    x = pc_loadings_filtered[top_genes_pc2, "PC1"],
                    y = pc_loadings_filtered[top_genes_pc2, "PC2"],
                    gene_symbol = top_genes_symbols_pc2
                  ),
                  aes(x = x, y = y, label = gene_symbol),
                  vjust = 1.5,
                  hjust = -0.5,
                  color = "darkorange",
                  size  = 2
                )
              # Customize axis labels to be bold and centered
              biplot <- biplot +
                theme(axis.title.x = element_text(face = "bold", hjust = 0.5),
                      axis.title.y = element_text(face = "bold", hjust = 0.5),
                      plot.title = element_text(face = "bold", hjust = 0.5))  
              # Save biplot
              ggsave("biplot.png", plot = biplot + 
                       theme(panel.background = element_rect(fill = "white")), 
                     width = 15, height = 8, units = "in", dpi = 300)
              
              library(plotly)
              
              # Plot each dataset with its corresponding color
              # Create separate datasets for each category
              upregulated <- top[top$logFC > input$logFC & top$adj.P.Val < input$adjPval, ]
              downregulated <- top[top$logFC < -input$logFC & top$adj.P.Val < input$adjPval, ]
              not_significant <- top[!(top$logFC > input$logFC & top$adj.P.Val < input$adjPval) & !(top$logFC < -input$logFC & top$adj.P.Val < input$adjPval), ]
              max_y <- max(-log10(c(upregulated$adj.P.Val, downregulated$adj.P.Val, not_significant$adj.P.Val)))
              
              # Create custom hover text including gene symbol, log2FC, and adjusted p-value
              up_text <- paste("Gene Symbol: ", upregulated$gene_symbol, "<br>",
                               "log2FC: ", format(round(upregulated$logFC, 2), nsmall = 2), "<br>",
                               "Adjusted P-value: ", format(round(upregulated$adj.P.Val, 3), nsmall = 3))
              
              down_text <- paste("Gene Symbol: ", downregulated$gene_symbol, "<br>",
                                 "log2FC: ", format(round(downregulated$logFC, 2), nsmall = 2), "<br>",
                                 "Adjusted P-value: ", format(round(downregulated$adj.P.Val, 3), nsmall = 3))
              
              not_significant_text <- paste("Gene Symbol: ", not_significant$gene_symbol, "<br>",
                                            "log2FC: ", format(round(not_significant$logFC, 2), nsmall = 2), "<br>",
                                            "Adjusted P-value: ", format(round(not_significant$adj.P.Val, 3), nsmall = 3))
              
              # Create the plot with custom hover text
              plot <- plot_ly() %>%
                add_trace(data = upregulated, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                          marker = list(size = 4, color = 'red'), name = 'Upregulated', text = up_text, hoverinfo = "text") %>%
                add_trace(data = downregulated, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                          marker = list(size = 4, color = 'blue'), name = 'Downregulated', text = down_text, hoverinfo = "text") %>%
                add_trace(data = not_significant, x = ~logFC, y = ~-log10(adj.P.Val), type = 'scatter', mode = 'markers',
                          marker = list(size = 2, color = 'black'), name = 'Not significant', text = not_significant_text, hoverinfo = "text") %>%
                layout(
                  title = list(text = paste0("Interactive Volcano Plot"), y = 0.97),  
                  xaxis = list(title = "log2 Fold Change", showgrid = FALSE, zeroline = FALSE), 
                  yaxis = list(title = "-log10(Adjusted P-value)", showgrid = FALSE, linecolor= "black"), 
                  hovermode = "closest",
                  legend = list(
                    orientation = "v",
                    x = 0.9,
                    y = 0.2,
                    xanchor = "left",
                    yanchor = "top",
                    bgcolor = "rgba(0,0,0,0)",
                    bordercolor = "rgba(0,0,0,0)",  
                    borderwidth = 2
                  ),
                  margin = list(l = 50, r = 50, b = 50, t = 50),  
                  font = list(family = "Arial", size = 14, color = "black", weight = "bold"),
                  showlegend = TRUE,
                  yaxis2 = list(range = c(0, max_y + 1)), 
                  print("volcano is ready")
                )
              
              htmlwidgets::saveWidget(plot, file = "volcanoplot.html")
              
              # Heatmap
              
              library(RColorBrewer)
              library(heatmaply)
              
              sel.top<- subset(top, top$adj.P.Val<input$adjPval & abs(top$logFC)>input$logFC)
              sel.top<- annotf(sel.top, y)
              sel.top<- unFilt(sel.top$gene_symbol, sel.top)
              sel.top<-sel.top[is.na(sel.top$gene_symbol)==FALSE,]
              sel.top<- sel.top[order(sel.top$logFC, decreasing = TRUE),]
              f2_e.sel<- f2_e[rownames(sel.top), ]

              keys<-sel.top$probeID
              values<-sel.top$gene_symbol
              l<- list()
              for(i in 1:length(keys)){
                l[keys[i]]<- values[i]
              }
              rownames(f2_e.sel)<- unlist(l[rownames(f2_e.sel)], use.names = F)
              yellowblackblue <- colorRampPalette(c("dodgerblue", "black", "gold"))(n = 100)
              data<- matrix(data = f2_e.sel, nrow = nrow(f2_e.sel))
              hover_text <- matrix(paste("AveExpr:", round(data, 4)), 
                                   nrow = nrow(data), ncol = ncol(data))
              
              heatdata<- f2_e.sel
              colnames(heatdata)<-empty_up$sample
              
              heat<- heatmaply(heatdata, distfun = "pearson", hclust_method = NA, 
                               colors = yellowblackblue, row_dend_left = T,
                               margins = c(50,50,70,0),
                               main = "Heatmap", scale = "row", Colv = F,
                               label_names = c("Gene", "Sample", "Value"),
                               custom_hovertext = hover_text,
                               file = "./heatmap.html")
              #browseURL("./heatmap.html")
              
              files_to_zip<- c("volcanoplot.html",
                               "boxplot.png",
                               "histogram.png",
                               "scree_plot.png",
                               "pca_scores_plot.png",
                               "biplot.png",
                               "heatmap.html"
              )
              zip_filename<- "plots.zip"
              zip(zip_filename, files_to_zip)
              
              output$down_plots <- downloadHandler(
                filename = function() {
                  zip_filename
                },
                content = function(file) {
                  path_to_zip <- zip_filename
                  file.copy(path_to_zip, file)
                }
              )
              
              ###Filtering for NA geneNames
              sign<- sign[is.na(sign$gene_symbol)==FALSE,]
              print(paste0("nrow(sign) after removing NAs is ",nrow(sign)))
              
              ###Filtering for unique geneNames
              sign_un<- unFilt(sign$gene_symbol, sign)
              print(paste0("DEGs table after unFilt:",dim(sign_un)))
              if (nrow(sign_un)==0){
                ###Message for zero DE genes
                print("nrow of sign_un is 0 - no DEGs")
                shinyjs::hideElement ("wait_msg3")
                shinyjs::hideElement ("wait_msg2")
                shinyjs::showElement ("msg")
                shinyjs::showElement ("noDEGs")
              }
              else {
                ###Rendering and downloading DE genes list
                print(paste0("If the number of rows of sign_un is >=1: ", nrow(sign_un)))
                sign_un <- sign_un[order(sign_un$adj.P.Val),]
                degs <- sign_un
                degs_t<- sign_un
                sign_un$gene_symbol <- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/gene?term=", sign_un$gene_symbol,
                                                                  "[Gene Name] AND ", as.character(unique(org)), 
                                                                  "[Organism]"),
                                              "' target='_blank'>", sign_un$gene_symbol,"</a>")
                shinyjs::showElement ("DEGs")
                sign_un<- DT::datatable( data = sign_un,                       
                                         editable = FALSE, 
                                         colnames = c("probeID","gene symbol","logFC","AveExpr",
                                                      "t","P.Value","adj.P.Val","B"),
                                         rownames = FALSE,
                                         escape = FALSE,
                                         extensions = c('Scroller' , 'Buttons'),
                                         options = list (
                                           autoWidth = FALSE,
                                           dom = 'Blftrip',
                                           deferRender = T,
                                           scrollY = 750,
                                           scroller = T,
                                           buttons = list(list(extend = 'colvis', columns = c(2, 3, 4, 5, 6, 7)),
                                                          "copy", "print")
                                         )
                )%>% formatSignif(c('logFC','AveExpr','t','P.Value','adj.P.Val','B'),3)
                
                gene_link <- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/gene?term=", degs$gene_symbol,
                                                        "[Gene Name] AND ", as.character(unique(org)), 
                                                        "[Organism]"),
                                    "' target='_blank'>", degs$gene_symbol,"</a>")
                degs$gene_symbol <- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/gene?term=", degs$gene_symbol,
                                                               "[Gene Name] AND ", as.character(unique(org)), 
                                                               "[Organism]"),
                                           "' target='_blank'>", degs$gene_symbol,"</a>")
                degs <- cbind(degs, url = gene_link)
                
                output$DEGs <- DT::renderDataTable(sign_un, server = FALSE)
                output$down_csv <- downloadHandler(
                  filename = function (){ paste0("GSE",input$GSE,"_DE_genes.csv")},
                  content = function (temp) { write.csv(degs, temp, row.names = FALSE)})
                output$down_tsv <- downloadHandler(
                  filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                  content = function (temp) { write.table(degs, temp, quote = FALSE, sep = "\t",row.names = FALSE) })
                output$down_pdf <- downloadHandler( 
                  filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                  content = function ( temp ) {
                    createPdf(degs)
                    file.copy ("test.pdf", temp)
                  },
                  contentType = "application/pdf")
                
                shinyjs::hideElement ("wait_msg2")
                shinyjs::showElement ("msg")
                shinyjs::hideElement ("wait_msg3")
                shinyjs::showElement ("download")
                shinyjs::enable ("download")
                shinyjs::enable ("changes_done")
                shinyjs::enable ("submit_comp")
                shinyjs::enable ("anls")
                shinyjs::hideElement ("wait_msg5")
                shinyjs::showElement ("WebGestaltR_but")

                library(WebGestaltR)
                
                observeEvent(input$WebGestaltR,{
                  shinyjs::showElement("WebGestaltR_but")
                  shinyjs::disable("WebGestaltR_but")
                  shinyjs::showElement("choosing_organism")
                  output$chooseOrganism<- renderUI({
                    selectInput( inputId= "chooseOrganism",
                                 label = "Choose the organism for the analysis",
                                 choices = listOrganism()
                    )
                  })
                  shinyjs::enable("submit_Organism")
                  print("Choosing organism")
                  
                  observeEvent(input$submit_Organism,{
                    shinyjs::disable("submit_Organism")
                    shinyjs::showElement("choosing_RefSet")
                    organismWebGestaltR<- input$chooseOrganism
                    output$chooseRefSet<- renderUI({
                      selectInput( inputId = "chooseRefSet",
                                   label = "Choose the reference set for the analysis",
                                   choices = listReferenceSet(organism = organismWebGestaltR))})
                    print("Choosing RefSet")
                    
                    observeEvent(input$submit_RefSet,{
                      shinyjs::disable("submit_RefSet")
                      shinyjs::hideElement("wait_msg5")
                      shinyjs::showElement("wait_msg6")
                      shinyjs::hideElement("WebGestalt_down")
                      shinyjs::hideElement("WebGestalt_res")
                      shinyjs::hideElement("WebGestaltR_but")
                      shinyjs::hideElement("choosing_organism")
                      shinyjs::hideElement("choosing_RefSet")

                      library(zip)

                      print("Running WebGestaltR")
                      genesymbol <- degs_t$gene_symbol
                      genesymbol <- as.character(genesymbol)
                      projectName <- c("WebGestalt")
                      outputDir <- currentDir
                      WebGestaltR(enrichMethod = "ORA", organism = input$chooseOrganism,
                                  enrichDatabase = c("geneontology_Biological_Process", "geneontology_Molecular_Function",
                                                     "geneontology_Cellular_Component"),
                                  interestGene= genesymbol, interestGeneType="genesymbol", 
                                  collapseMethod="mean", referenceSet= input$chooseRefSet,
                                  minNum=10, maxNum=500, sigMethod="fdr", fdrMethod="BH", fdrThr=0.05, 
                                  reportNum=20, isOutput= T, outputDirectory=outputDir,
                                  projectName=projectName, hostName= "https://www.webgestalt.org")
                      HTMLfile<- list.files(path = paste0(outputDir,"/Project_",projectName), pattern = ".html")
                      print(paste0("HTML file is ",HTMLfile))
                      projectPath<- paste0(outputDir,"/","Project_", projectName)
                      print(paste0("project path is ", projectPath))
                      file.copy(from = paste0(projectPath,"/",HTMLfile),to = paste0(outputDir,"/www"))
                      shinyjs::hideElement("wait_msg6")
                      shinyjs::showElement("WebGestalt_down")
                      shinyjs::showElement("WebGestalt_res")
                      output$WebGestalt <- renderUI({tags$iframe(style="height:800px; width:100%", src=HTMLfile)})
                      Webgestalt_files<- list.files(path = paste0(outputDir,"/Project_", projectName), 
                                                    recursive = TRUE, full.names = TRUE)
                      zip::zipr(zipfile = "WebGestalt.zip", files = Webgestalt_files, include_directories=TRUE)
                      if (file.exists("WebGestalt.zip")) {print ("zip is ready!!!")}
                      # the zip file is in working directory
                      else {print("Problem while zipping")}
                    })
                    output$down_WebGestalt <- downloadHandler(
                      filename = function() {
                        paste0("WebGestalt_Results.zip","")
                      },
                      content = function(temp) {
                        file.copy(from = paste0(getwd(),"/WebGestalt.zip"), to = temp)  
                      })
                  }) #submitOrganism
                }) #WebGestlatR
              } #else sto Rendering
            } #!is.element
            else { #sto if(!is.element(y,set))
              shinyjs::showElement("no_annotation")
              
              ###Rendering and downloading DE genes list WITHOUT ANNOTATION 
              sign<- no_annotf(sign)
              dim(sign)
              sign2 <- sign[order(sign$adj.P.Val),]
              degs2<- sign2
              shinyjs::showElement ("DEGs")
              sign2<- DT::datatable( data = sign2,                       
                                     editable = FALSE, 
                                     colnames = c("probeID","logFC","AveExpr","t","P.Value","adj.P.Val","B"),
                                     rownames = FALSE,
                                     extensions = c('Scroller' , 'Buttons'),
                                     options = list (
                                       autoWidth = FALSE,
                                       dom = 'Blftrip',
                                       deferRender = T,
                                       scrollY = 750,
                                       scroller = T,
                                       buttons = list(list(extend = 'colvis', columns = c(1, 2, 3, 4, 5, 6)),
                                                      "copy", "print")
                                     )
              )%>% formatSignif(c('logFC','AveExpr','t','P.Value','adj.P.Val','B'),3)
              
              output$DEGs <- DT::renderDataTable(sign2, server = FALSE)
              output$down_csv <- downloadHandler(
                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.csv")},
                content = function (temp) { write.csv(degs2, temp, row.names = FALSE)})
              output$down_tsv <- downloadHandler(
                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                content = function (temp) { write.table(degs2, temp, quote = FALSE, sep = "\t",row.names = FALSE) })
              output$down_pdf <- downloadHandler( 
                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                content = function ( temp ) {
                  createPdf(degs2)
                  file.copy ("test.pdf", temp)
                },
                contentType = "application/pdf")
              
              shinyjs::hideElement ("wait_msg2")
              shinyjs::showElement ("msg")
              shinyjs::hideElement ("wait_msg3")
              shinyjs::showElement ("download")
              shinyjs::enable ("download")
              shinyjs::enable ("changes_done")
              shinyjs::enable ("submit_comp")
              shinyjs::enable ("anls")
              shinyjs::hideElement("wait_msg5")
              shinyjs::showElement("msg_no_annot")
            } #no annotation
          } #else - annotate
          
          #Deleting files from working directory
          print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
          unlink(file.path(getwd(),c("*.CEL", "*.tar", "*.gz","*.txt.gz","*.txt"))) 
          print(".CEL, *.tar, *.gz, *.txt.gz, and *.txt are gone")
          print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
          
        }) #anls
      }) #submit_comp
    }) #changes_done
  }) #up_go
}) #server
