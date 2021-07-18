#Loading packages
library(shiny)
library(shinyjs)
library(shinyBS)
library(DT)
library(RCurl)
library(BiocManager)
library(GEOquery)
library(limma)
library(oligo)
library(annotate)

source('helpers.R')

#setwd('../GSE')

# Define server logic required to download a GSE matrix
shinyServer(function(input, output) {
    shinyjs:: useShinyjs()  # Include shinyjs
    shinyjs::showElement("wait_msg5")
    #Deleting files from working directory
    print(paste0("working directory is ", getwd()))
    #GSEdir<- getwd()
    #print(paste0("GSEdir is ", GSEdir))
    #unlink(paste0(GSEdir,"*"),recursive = TRUE)
    currentDir<- "./"
    ProjectDir<- dir (path = currentDir, pattern ="^Project", full.names = TRUE)
    print(paste0("ProjectDir is ", ProjectDir))
    print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
    #unlink(paste0(ProjectDir,"*"), recursive = TRUE)
    unlink(file.path(getwd(),c("*.CEL", "*.gz", "*.txt.gz","*.txt", "*.tar"))) 
    print("*.CEL, *.tar, *.gz, *.txt.gz and *.txt are gone")
    print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
    wwwDir<- "./www"    
    unlink(file.path(wwwDir,c ("*.html")))

    #changing the file size limit
    options(shiny.maxRequestSize = 25*1024^2) #25MB per file is allowed
    
    #preventing the shiny app from being grayed out
    autoInvalidate <- reactiveTimer(10000)
    observe({
        autoInvalidate()
        cat(".")
    })
    
    #View user's guide
    output$users_guide_msg <- renderUI({tagList(tags$p(tags$b(tags$h2("User's guide"))))})
    output$users_guide_view <- renderUI({
        tags$iframe(style="height:800px; width:100%", src="Users_guide.pdf")
    })
    
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
        
        mandatoryFilled <- vapply(fieldsMandatory, function(x) {
            !is.null(input[[x]]) && input[[x]] != ""
        },
        logical(1))
        mandatoryFilled <- all(mandatoryFilled)
        
        mandatory_upload_Filled <- vapply(fieldsMandatory_up, function(x) {
            !is.null(input[[x]])},
            logical(1))
        mandatory_upload_Filled <- all(mandatory_upload_Filled)
        
        #enable/disable the submit button
        shinyjs::toggleState(id = "go", condition = mandatoryFilled)
        shinyjs::toggleState(id = "up_go", condition = mandatory_upload_Filled)
    })
    
    #after submission
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
        
        #shinyjs::showElement ("choosing_platform")
        
        
        #GSE accession number validation and downloading
        accession<- input$GSE
        print(paste0("The length of input$GSE is ",nchar(accession)))
        ftp<- createFtp(accession)
        
        if (validation(ftp)) {
            link <-createLink(input$GSE)
            gse<- createGSE2(input$GSE)
            output$url <- renderUI({tagList(tags$p("URL link:", link))})

            shinyjs::showElement ("goToNextTab")
            #Choosing platform
            plat_choices <- c()
            choose_platform <- FALSE
            shinyjs::showElement("choosing_platform")
            shinyjs::hideElement ("wait_msg4")
            if (length(gse)>0) {
                choose_platform <- TRUE
                for (i in 1:length(gse)) {
                    plat_choices <- c(plat_choices, annotation(gse[[i]])) #class(plat_choices) >>> character
                }
            }
            plat_df <- as.data.frame(plat_choices)
            for (i in 1:length(plat_choices)) {
                row.names(plat_df)[i]<- plat_choices[i]
            }
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
            
            shinyjs::disable ("discard")
            
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
            } else {
                for (i in 1:(ncol(pdata))) {
                    for (j in 1:(nrow(pdata))) {
                        cont <- as.character(strsplit(paste(pdata[j, i]), ":")[[1]])
                        newpdata[j, i] <- cont[length(cont)]
                    }
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
                               choices = (colnames(newpdata))
                               )
                })
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
                        if ("replicate" %in% colnames(newpdata)) {
                            
                        }
                        else {
                            replicate<- matrix(data = integer(), nrow = dim(newpdata)[1], ncol = length(unique(newpdata[, index])))
                            row.names(replicate) <- row.names(newpdata)
                            colnames(replicate) <- unique(newpdata[, index])
                            replicate[1,1]<-1
                            for (i in 2:dim(newpdata)[1]){
                                c<- which(colnames(replicate)==newpdata[i, index], arr.ind = T)
                                if (! is.na(replicate[i-1,c])) {
                                    replicate[i,c]<- replicate[i-1,c]+1
                                }
                                else {
                                    replicate[i,c]<-1
                                }
                            }
                            replicate
                            replicate<- as.integer(replicate)
                            replicate<- replicate[!is.na(replicate)]
                            newpdata <- cbind(newpdata, replicate= replicate)
                        }
                        ###Creating targets.txt
                        write.table(newpdata, file = "targets.txt", row.names = FALSE, quote = FALSE, sep = "\t")
                        if (file.exists("targets.txt")) {
                            print("targets.txt is ready!!!")
                        }
                        
                        #Back to top button
                        shinyjs::enable("toTop")
                        observeEvent(input$toTop, {
                            js$toTop();
                        })
                        
                        #run the analysis
                    observeEvent(input$anls,{
                        shinyjs::disable ("toTop")
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
                        targets<- readTargets(file="targets.txt", sep="\t")
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
                        pl_ctr <- length(which(newpdata[,index] %in% input$ctr)) ###plh8os control samples
                        pl_tr <- length(which(newpdata[,index] %in% input$tr)) ###plh8os treated samples
                        p<- (pl_tr/(pl_tr+pl_ctr))
                        print(paste("pl_ctr is", pl_ctr, sep = " "))
                        print(paste("pl_tr is", pl_tr, sep = " "))
                        print(paste("p is", p, sep = " "))
                        f2_Data_rma<-filter1(Data_rma,p,log2(100))  ###p=plh8os treated/synoliko plh8os samples
                        dim(f2_Data_rma)
                        f2_e<- filter2(e, p, log2(100))
                        dim(f2_e)
                        ###NonSpecific Filtering-shorth
                        row.mean <- esApply(Data_rma,1,mean) 
                        sh <- shorth(row.mean)
                        hist(e)
                        abline(v=sh, col="red")
                        f3_Data_rma <- Data_rma[row.mean >=sh,]
                        dim(f3_Data_rma)
                        ###Filtering Results - Comparison
                        v<- c(dim(Data_rma)[1], dim(f2_Data_rma)[1], dim(f3_Data_rma)[1])
                        filt.matrix<- data.frame(features=v, filt.method = c("before filt", "pOverA", "shorth"))  
                        print(filt.matrix)
                        ###Quality Assessmament- Checking for Batch Effect
                        plotMDS(e, labels = targets$treatment, col=as.numeric(as.factor(targets$replicate)))
                        plotMDS(f2_e, labels = targets$treatment, col=as.numeric(as.factor(targets$replicate)))
                        ###Linear Modelling
                        cond<-pdataRel[,index]
                        batch<-as.factor(targets$replicate)
                        design<-model.matrix(~0 + cond + batch)
                        x_index <- grep(input$ctr, colnames(design), ignore.case = T) ###
                        y_index <- grep(input$tr, colnames(design), ignore.case = T) ###
                        colnames(design)<-gsub("cond","",colnames(design))
                        colnames(design)<-gsub("atch","",colnames(design))
                        colnames(design)<- make.names(colnames(design))
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
                        vennDiagram(res, include = "up", show.include = T)
                        vennDiagram(res, include = "down", show.include = T)
                        vennDiagram(res, include = "both", show.include = T)
                        
                        ###Keeping DE genes to a data.frame
                        sign<-top[which(abs(top$logFC) > input$logFC & top$adj.P.Val < input$adjPval),]
                        dim(sign)

                        ###Annotation using annotate
                        dict <- read.table("./plat_dict_full.tsv",
                                           header = T, sep = "\t", colClasses = "character")
                        gse_acc <- paste0("GSE", input$GSE)
                        gse<- getGEO(gse_acc, destdir = getwd(), GSEMatrix = T)
                        a <- annotation(gse[[idx]])
                        p<- grep(a, dict$Accession, ignore.case = T, fixed = T)
                        x<- dict$db[p]
                        y<-x[1]
                        
                        set<- c("-", "?", "")
                        if(!is.element(y,set)){ #checking if annotation file exists
                            
                            if(!require(y, character.only = T)) {
                                BiocManager::install(y, character.only = T, dependencies = T, ask=FALSE)
                                library(y, character.only = T)
                            }
                            
                            if(dim(sign)==0) {
                                shinyjs::hideElement ("wait_msg3")
                                shinyjs::hideElement ("wait_msg2")
                                shinyjs::showElement ("msg")
                                shinyjs::showElement ("noDEGs")
                            }
                            else {
                                sign<- annotf1(sign, y)
                                dim(sign)
                                ###Filtering for NA geneNames
                                sign<- sign[is.na(sign$gene_symbol)==FALSE,]
                                dim(sign)
                                ###Filtering for unique geneNames
                                sign_un<- unFilt(sign$gene_symbol, sign)
                                print(paste0("DEGs table:",dim(sign_un)))
                                
                                ###Rendering and downloading DE genes list
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
                                degs$gene_symbol <- gene_link#paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", degs$gene_symbol),
                                                    #       "' target='_blank'>", degs$gene_symbol,"</a>")
                                degs <- cbind(degs, url = gene_link)
                                
                                output$DEGs <- DT::renderDataTable(sign_un, server = FALSE)
                                
                                #output$down_degs_t <- downloadHandler(
                                #    filename = function(){paste0("GSE", input$GSE,"_DE_genes.csv")},
                                #    content = function(temp) {write.csv(degs_t, temp, row.names = FALSE)}
                                #)

                                output$down_csv <- downloadHandler(
                                    filename = function (){ paste0("GSE",input$GSE,"_DE_genes.csv")},
                                    content = function (temp) { write.csv(degs, temp, row.names = FALSE)}
                                )
                                output$down_tsv <- downloadHandler(
                                    filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                                    content = function (temp) { write.table(degs, temp, quote = FALSE, sep = "\t",row.names = FALSE) }
                                )
                                output$down_pdf <- downloadHandler( 
                                    filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                                    content = function ( temp ) {
                                        createPdf(degs)
                                        file.copy ("test.pdf", temp)
                                    },
                                    contentType = "application/pdf"
                                )
                                
                                shinyjs::hideElement ("wait_msg2")
                                shinyjs::showElement ("msg")
                                shinyjs::hideElement ("wait_msg3")
                                shinyjs::showElement ("download")
                                shinyjs::enable ("download")
                                shinyjs::enable ("changes_done")
                                shinyjs::enable ("submit_comp")
                                shinyjs::enable ("anls")
                                shinyjs::disable ("toTop")
                                shinyjs::hideElement ("wait_msg5")
                                shinyjs::showElement ("WebGestaltR_but")
                                
                                ###WebGestaltR
                                #if(!require(WebGestaltR)) {
                                #    install.packages('WebGestaltR', dependencies = TRUE)
                                    library(WebGestaltR)
                                #}
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
                                                     choices = listReferenceSet(organism = organismWebGestaltR)
                                                     )
                                            })
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
                                            #if(!require(zip)) {
                                            #    install.packages('zip', dependencies = TRUE)
                                                library(zip)
                                            #}
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
                                                        projectName=projectName, hostName= "http://www.webgestalt.org"
                                            )
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
                                            if (file.exists("WebGestalt.zip")) {
                                                print ("zip is ready!!!")
                                                # the zip file is in working directory, here C:/Users/Anna/Desktop/PhD/DExplore/GSE
                                            }
                                            else {print("Problem while zipping")}
                                        })
                                        output$down_WebGestalt <- downloadHandler(
                                            filename = function() {
                                                paste0("WebGestalt_Results.zip","")
                                            },
                                            content = function(temp) {
                                                file.copy(from = paste0(getwd(),"/WebGestalt.zip"), to = temp)  
                                            }
                                        )
                                        })
                                    })
                            }
                        }
                        else {
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
                                content = function (temp) { write.csv(degs2, temp, row.names = FALSE)}
                            )
                            output$down_tsv <- downloadHandler(
                                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                                content = function (temp) { write.table(degs2, temp, quote = FALSE, sep = "\t",row.names = FALSE) }
                            )
                            output$down_pdf <- downloadHandler( 
                                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                                content = function ( temp ) {
                                    createPdf(degs2)
                                    file.copy ("test.pdf", temp)
                                },
                                contentType = "application/pdf"
                            )
                            
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
                        }
                        #Deleting files from working directory
                        #unlink(paste0(GSEdir,"*"), recursive = TRUE)
                        print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
                        unlink(file.path(getwd(),c("*.CEL", "*.tar", "*.gz", "*.txt.gz","*.txt"))) 
                        print("*.CEL, *.tar, *.gz, *.txt.gz and *.txt are gone")
                        print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
                        #Back to top button
                        shinyjs::enable("toTop")
                        observeEvent(input$toTop, {
                            js$toTop();
                        })
                        
                    })
                    })
                    })
                })
            })
        
    
        }
        
        else {
            output$falseGSE <- renderText({
                paste("GSE",input$GSE, " is not a valid accesion number. Please, check again.", sep = "")}) 
            shinyjs::showElement ("falseGSE")
            print("validation==FALSE")
        }
        shinyjs::hideElement ("submit_msg")
        shinyjs::hideElement ("wait_msg")
    }
    )
    
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
        shinyjs::disable ("discard")
        
        
        #Editing Data Description form
        observeEvent(input$df_up_cell_edit,{
            info <- input$df_up_cell_edit
            i <- info$row
            j <- info$col+1
            v <- info$value
            empty_up[i,j]<<- DT::coerceValue(v,empty_up[i,j])
            shinyjs::enable ("discard")
        })
        
        #Saving changes - Creating targets.txt
        observeEvent(input$changes_done,{
            shinyjs::disable("changes_done")
            print("changes_done")
            ###comparisons
            if(length(unique(empty_up$treatment))>1)     { a<- empty_up$treatment }      else { a <-NULL }
            if(length(unique(empty_up$duration))>1)      { b<- empty_up$duration }       else { b <-NULL }
            if(length(unique(empty_up$concentration))>1) { c<- empty_up$concentration }  else { c <-NULL }
            abc <- paste(a, b, c, sep = "-")
            if ( (length(unique(abc))>=2) & (length(unique(empty_up$replicate))>=2) ){
                df_comp <- create_df_comp(unique(abc))
                print(df_comp)
                m_comp <- create_m_comp(df_comp)
                print(m_comp)
                shinyjs::hideElement ("no_comparison")
                shinyjs::showElement ("comparison")
                output$compBut <- renderUI ({ 
                    radioButtons ( inputId = "comp2",
                                   label = "Select the comparison",
                                   choices = (m_comp)
                    ) 
                })
                shinyjs::enable("submit_comp")
                
            }
            else  {
                shinyjs::showElement ("no_comparison")
                output$comp <- renderText( c ("Please, fill in the form correctly.",
                                              "Columns 'treatment' and 'replicate' should be filled out!",
                                              "You may use information from GEO's webpage."))
                
            }
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
            write.table(empty_up, file = "targets.txt", row.names = FALSE, quote = FALSE, sep = "\t")
            
            observeEvent(input$submit_comp, {
                print("comp_submitted")
                print(input$comp2)
                print(class(input$comp2))
                index <- which(m_comp==input$comp2)
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
                    
                    ###Reading targets.txt
                    targets<- readTargets(file="targets.txt", sep="\t")
                    
                    ###Reading Affy
                    Data<- read.celfiles(input$uploadFiles$datapath)
                    Data_rma<-rma(Data)
                    dim(Data_rma)
                    e <- exprs(Data_rma)
                    print(dim(e))
                    ###NonSpecific Filtering-pOverA
                    pl_tr<-length(subset(abc, abc==df_comp[index,3])) ###plh8os treated samples
                    pl_ctrl<-length(subset(abc, abc==df_comp[index,1])) ###plh8os control samples
                    p<- (pl_tr/(pl_tr+pl_ctrl))
                    print(paste0("p is ",p))
                    f2_Data_rma<-filter1(Data_rma,p,log2(100))  ###p=plh8os treated/plh8os control samples
                    dim(f2_Data_rma)
                    f2_e<- filter2(e, p, log2(100))
                    dim(f2_e)
                    ###NonSpecific Filtering-shorth
                    row.mean <- esApply(Data_rma,1,mean) 
                    sh <- shorth(row.mean)
                    hist(e)
                    abline(v=sh, col="red")
                    f3_Data_rma <- Data_rma[row.mean >=sh,]
                    dim(f3_Data_rma)
                    ###Filtering Results - Comparison
                    v<- c(dim(Data_rma)[1], dim(f2_Data_rma)[1], dim(f3_Data_rma)[1])
                    filt.matrix<- data.frame(features=v, filt.method = c("before filt", "pOverA", "shorth"))  
                    print(filt.matrix)
                    ###Quality Assessment- Checking for Batch Effect
                    plotMDS(e, labels = targets$treatment, col=as.numeric(as.factor(targets$replicate)))
                    plotMDS(f2_e, labels = targets$treatment, col=as.numeric(as.factor(targets$replicate)))
                    ###Linear Modelling
                    abc<-gsub("-","_",abc)
                    cond<-as.factor(abc)
                    batch<-as.factor(targets$replicate)
                    design<-model.matrix(~0 + cond + batch)
                    colnames(design)<-gsub("cond","",colnames(design))
                    colnames(design)<-gsub("atch","",colnames(design))
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
                    vennDiagram(res, include = "up", show.include = T)
                    vennDiagram(res, include = "down", show.include = T)
                    vennDiagram(res, include = "both", show.include = T)
                    
                    ###Keeping DE genes to a data.frame
                    sign<-top[which(abs(top$logFC) > input$logFC & top$adj.P.Val < input$adjPval),]
                    dim(sign)
                    
                    ###Annotation me annotate
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
                        
                        if(dim(sign)==0) {
                            shinyjs::hideElement ("wait_msg3")
                            shinyjs::hideElement ("wait_msg2")
                            shinyjs::showElement ("msg")
                            shinyjs::showElement ("noDEGs")
                        }
                        else {
                            sign<- annotf1(sign, y)
                            dim(sign)
                            ###Filtering for NA geneNames
                            sign<- sign[is.na(sign$gene_symbol)==FALSE,]
                            dim(sign)
                            ###Filtering for unique geneNames
                            sign_un<- unFilt(sign$gene_symbol, sign)
                            print(paste0("DEGs table:",dim(sign_un)))
                            
                            ###Rendering and downloading DE genes list
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
                                content = function (temp) { write.csv(degs, temp, row.names = FALSE)}
                            )
                            output$down_tsv <- downloadHandler(
                                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                                content = function (temp) { write.table(degs, temp, quote = FALSE, sep = "\t",row.names = FALSE) }
                            )
                            output$down_pdf <- downloadHandler( 
                                filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                                content = function ( temp ) {
                                    createPdf(degs)
                                    file.copy ("test.pdf", temp)
                                },
                                contentType = "application/pdf"
                            )
                            
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
                            ###WebGestaltR
                            #if(!require(WebGestaltR)) {
                            #    install.packages('WebGestaltR', dependencies = TRUE)
                                library(WebGestaltR)
                            #}
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
                                                     choices = listReferenceSet(organism = organismWebGestaltR)
                                        )
                                    })
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
                                        #if(!require(zip)) {
                                        #    install.packages('zip', dependencies = TRUE)
                                            library(zip)
                                        #}
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
                                                    projectName=projectName, hostName= "http://www.webgestalt.org"
                                        )
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
                                        if (file.exists("WebGestalt.zip")) {
                                            print ("zip is ready!!!")
                                            # the zip file is in working directory, here C:/Users/Anna/Desktop/PhD/DExplore/GSE
                                        }
                                        else {print("Problem while zipping")}
                                    })
                                    output$down_WebGestalt <- downloadHandler(
                                        filename = function() {
                                            paste0("WebGestalt_Results.zip","")
                                        },
                                        content = function(temp) {
                                            file.copy(from = paste0(getwd(),"/WebGestalt.zip"), to = temp)  
                                        }
                                    )
                                })
                            })
                        }
                    }
                    #shinyjs::hideElement("wait_msg5")
                    #shinyjs::showElement("wait_msg6")
                    #shinyjs::hideElement("WebGestalt_down")
                    #shinyjs::hideElement("WebGestalt_res")
                    else {
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
                            content = function (temp) { write.csv(degs2, temp, row.names = FALSE)}
                        )
                        output$down_tsv <- downloadHandler(
                            filename = function (){ paste0("GSE",input$GSE,"_DE_genes.tsv")},
                            content = function (temp) { write.table(degs2, temp, quote = FALSE, sep = "\t",row.names = FALSE) }
                        )
                        output$down_pdf <- downloadHandler( 
                            filename = function (){ paste0("GSE",input$GSE,"_DE_genes.pdf")},
                            content = function ( temp ) {
                                createPdf(degs2)
                                file.copy ("test.pdf", temp)
                            },
                            contentType = "application/pdf"
                        )
                        
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
                    }
                    
                    #Deleting files from working directory
                    #unlink(paste0(GSEdir,"*"), recursive = TRUE)
                    print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
                    unlink(file.path(getwd(),c("*.CEL", "*.tar", "*.gz","*.txt.gz","*.txt"))) 
                    print(".CEL, *.tar, *.gz, *.txt.gz and *.txt are gone")
                    print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
                })
            })
        })
    })
    
    #resetting the form - not working
    #observeEvent(input$reset, {
        #print("reset")
        #shinyjs::reset ("form")
        #shinyjs::enable ("GSE")
        #shinyjs::showElement ("firstInput1")
        #shinyjs::showElement ("firstInput2")
        #shinyjs::hideElement ("firstRun")
        #shinyjs::hideElement ("wait_msg")
        #shinyjs::hideElement ("submit_msg")
        #shinyjs::hideElement ("goToNextTab")
        #shinyjs::hideElement ("falseGSE")
        #shinyjs::hideElement ("data_table")
        #shinyjs::hideElement ("comparison")
        #shinyjs::hideElement ("no_comparison")
        #shinyjs::hideElement ("changes")
        #shinyjs::hideElement ("btns")
        #shinyjs::hideElement ("run")
        #shinyjs::hideElement ("msg")
        #shinyjs::hideElement ("wait_msg2")
        #shinyjs::hideElement ("wait_msg3")
        #shinyjs::hideElement ("DEGs")
        #shinyjs::hideElement ("download")
        #shinyjs::hideElement ("noDEGs")
        #shinyjs::js$refresh()
        
        #Deleting files from working directory
        #unlink(paste0(GSEdir,"*"),recursive = TRUE)
        #print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
        #unlink(file.path(getwd(),c("*.CEL", "*.tar", "*.gz", "*.txt.gz","*.txt"))) 
        ###
        #unlink(file.path(getwd(),"*.CEL"))
        #print("*.CEL, *.tar, *.gz, *.txt.gz and *.txt are gone")
        #print(list.files(path = getwd(), all.files = T, full.names = T, include.dirs = T))
    #})
})