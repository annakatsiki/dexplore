#Global variables
fieldsMandatory <- c("GSE")
fieldsMandatory_up <- c("uploadFiles")
http <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE"

#Global functions
createFtp <- function(accession) {
  if (nchar(accession)<=3) {
    ftp <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE"
    ftp <- paste(ftp, accession, "/", sep="")
    print(paste0("FTP is ", ftp))
  }
  else if (nchar(accession)==4) {
    nnn <- substr(accession,1,1)
    ftp <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE"
    ftp <- paste(ftp, nnn, "nnn/GSE", accession, "/", sep="")
    print(paste0("FTP is ", ftp))
  }
  else if (nchar(accession)==5) {
    nnn <- substr(accession,1,2)
    ftp <- "ftp://nlm.nih.gov/geo/series/GSE"
    ftp <- paste(ftp, nnn, "nnn/GSE", accession, "/", sep="")
    print(paste0("FTP is ", ftp))
  }
  else if (nchar(accession)==6) {
    nnn <- substr(accession,1,3)
    ftp <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE"
    ftp <- paste(ftp, nnn, "nnn/GSE", accession, "/", sep="")
    print(paste0("FTP is ", ftp))
  }
  return(ftp)
}

validation <- function(ftp) {
  re<- GET (ftp)
  if (!is.null(re)) { 
    print("url exists")
    print("validation==TRUE")
    return(TRUE)
  }
  if (is.null(re)) { 
    print("url does not exist")
    return(FALSE)
  }
}

createLink <- function(dts) {
  https <- paste0(http, dts)
  link <- a (paste("GSE", dts, sep = ""), href=https, target="_blank")
  print(c("here is the link",paste(link)))
  return(link)
  print("link created")
  print(link)
}

createGSE<- function(dts) {
  accession<- as.character(paste0("GSE", dts))
  gse <- getGEO(accession, GSEMatrix = T, getGPL = F) 
  print("gse was created")
  return(gse)
}

create_uploaded_df <- function(fileNames){
  empty <- data.frame(sample= fileNames[,1],
                      treatment=paste("treatment"),
                      duration=paste("duration"),
                      concentration=paste("concentration"),
                      replicate=paste("replicate"),
                      stringsAsFactors = FALSE)
  print("empty table was created")
  return(empty)
}

create_uploaded_DataTable <- function(fileNames) {
  df <- DT::datatable(data.frame(sample= fileNames[,1],
                                 treatment=paste("treatment"),
                                 duration=paste("duration"),
                                 concentration=paste("concentration"),
                                 replicate=paste("replicate")),
                      editable = TRUE,
                      rownames = FALSE,
                      selection = "none",
                      extensions = c( 'Buttons', 'KeyTable', 'Responsive'),
                      options = list(autoWidth = FALSE,
                                     dom='Bfrtip',
                                     pageLength= 24,
                                     buttons = list(list(extend = 'colvis', columns = c(2, 3)))
                      )
  ) %>% formatStyle('sample',  color='black', backgroundColor = '#D1C4E9', fontWeight = 'bold')
  print("DataTable was created")
  return(df)
}

create_df_comp <- function(unique_table) {
  c <- numeric()
  i <- numeric()
  j <- 0
  df_comp <- data.frame(tr1=character(),
                        vs.=character(),
                        tr2=character(),
                        stringsAsFactors = FALSE)
  for (c in 1: length(unique_table)) {
    for (i in c:length(unique_table)) {
      if (c!=i){
        j <- j+1
        print(paste("c=",c,"i=",i, sep = " "))
        un_i <- unique_table[i]
        un_c <- unique_table[c]
        df_comp[j,1]<- un_c
        df_comp[j,2]<-"vs."
        df_comp[j,3]<- un_i
      }
    }
  }
  return(df_comp)
}

create_m_comp<- function(df_comp) {
  c<- numeric()
  m_comp<- matrix(ncol = 1)
  for (c in 1:nrow(df_comp)) {
    m_comp[c] <- paste( df_comp[c,1],df_comp[c,2],df_comp[c,3], sep = " ") 
  }
  return(m_comp)
}

filter1<- function(ExpressionSet, threshold1, threshold2) {
  if(require(genefilter)){
    set <- ExpressionSet
    f1 <- pOverA(threshold1,threshold2, na.rm = T)
    ffun <- filterfun(f1)
    filter <- genefilter(exprs(set),ffun)
    filtered <- set[filter,]
    return(filtered)
  }
}

filter2<- function(Matrix, threshold1, threshold2) {
  if(require(genefilter)){
    m <- Matrix
    f1 <- pOverA(threshold1,threshold2, na.rm = T)
    ffun <- filterfun(f1)
    filter <- genefilter(m,ffun)
    filtered <- m[filter,]
    return(filtered)
  }
}

annotf<- function(DataFrame, x) {  ###gia annotation me annotate
  if (require(x, character.only = T)) { 
    sign <- DataFrame
    ID <- row.names(sign)
    geneName <- data.frame(gene_symbol=getSYMBOL(ID,x))
    sign <- cbind("probeID"= row.names(sign), "gene_symbol"=geneName, sign)
    return(sign)
  }
}

no_annotf<- function(DataFrame) {  ###gia annotation me annotate
  sign <- DataFrame
  ID <- row.names(sign)
  sign <- cbind("probeID"= row.names(sign), sign)
  return(sign)
}

unFilt<- function(character, DataFrame) { ###filtering for unique genes, 1st argument: sign$gene_symbol
  c <- character
  sign <- DataFrame
  un <- unique(c)
  print(c("length(un)=",length(un)), quote = FALSE)
  y <- duplicated(c)
  sign_un <- sign[y==FALSE,]
  return(sign_un)
}

createPdf<- function(degs){
  maxrow <-  30
  idx <- seq(1,maxrow)
  npages <- ceiling(nrow(degs)/maxrow)
  pdf("test.pdf", height = 10, width = 12)
  grid.table(degs[idx,], rows = NULL)
  for(i in 2:npages) {
    grid.newpage();
    if(i*maxrow <= nrow(degs)) {
      idx <- seq ( 1+ ((i-1) * maxrow), i * maxrow)
    }
    else {
      idx <- seq ( 1+ ((i-1) * maxrow), nrow(degs))
    }
    grid.table( degs[idx, ], rows = NULL)
  }
  dev.off()
  print("pdf created")
}
