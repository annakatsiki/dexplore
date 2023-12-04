# dexplore
DExplore: an online tool for detecting differentially expressed genes from mRNA microarray experiments 
 
DExplore (www.dexplore.gr) is an online and easy-to-use web application for detecting differentially expressed genes using data from NCBI Gene Expression Omnibus (GEO), the international public repository that archives and freely distributes high-throughput functional genomics data submitted by the research community. 

DExplore can also be used locally to analyze data even if not submitted to NCBI’s GEO and allows privacy for users that have not yet submitted their data and wish to acquire a preliminary estimation of their experimental results. 

In addition, the application enables the user to proceed directly to functional enrichment analysis with the well-established online tool WebGestalt, using the results obtained from the preceding differential expression analysis of DExplore. 

DExplore is built using R programming language, shiny package and Bioconductor. 

A docker image can be found on https://hub.docker.com/r/akatsiki/dexplore
