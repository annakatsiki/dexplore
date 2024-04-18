#Loading libraries
library(shiny) 
library(shinyjs)
library(shinyBS)
library(RCurl)
library(DT)

source('helpers.R')

#Mandatory fields
labelMandatory <- function(label) {
  tagList(label,span("*", class = "mandatory_star")
  )}

appCSS <- ".mandatory_star { color: red; }"

# Define UI for application 
shinyUI(fluidPage ( 
  shinyjs::useShinyjs(),  # Include shinyjs
  shinyjs::inlineCSS(appCSS),
  tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
  style = "background-color:#B2DFDB",
  title= "DExplore - Differential Gene Expression Analysis",
  navbarPage( tags$b("DExplore - Differential Gene Expression Analysis"),
              selected = NULL,
              position = "fixed-top",
              inverse = TRUE,
              tabPanel ( title = "Data Input",
                         br(),
                         br(),
                         br(),
                         br(),
                         h3("You can use either an NCBI GEO's accession number OR your own .CEL files", style="text-align:center"),
                         splitLayout ( cellWidths = c ("80%", "20%"), 
                                       column ( 2,
                                                div ( id = "form",
                                                      style = "text-align:center",
                                                      textInput( inputId = "GSE", labelMandatory("Enter a valid GSE accession number")),
                                                      tags$head(tags$style(HTML("#go{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;text-align:center;}"))),
                                                      actionButton( inputId = "go", label = "Submit"),
                                                      p(span("*", class = "mandatory_star")," This field is mandatory")
                                                ),
                                                bsTooltip(id= "GSE", 
                                                          title = "Please, enter only the number (e.g. for GSE41827, enter 41827)",
                                                          placement = "bottom",
                                                          trigger = "hover")),
                                       column ( 10,
                                                div( id = "upload",
                                                     style = "text-align:center",
                                                     fileInput( inputId = "uploadFiles", labelMandatory("Upload your own .CEL files"),
                                                                multiple = TRUE, accept = ".CEL"),
                                                     tags$head(tags$style(HTML("#up_go{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;}"))),
                                                     actionButton( inputId = "up_go", label = "Submit"),
                                                     p(span("*", class = "mandatory_star")," This field is mandatory")
                                                ))),
                         hr(),
                         h4("Please, refresh DExplore between analyses", style="text-align:center"),
                         hr(),
                         
                         column(12,
                                shinyjs::hidden(
                                  div (id = "submit_msg",
                                       style = "text-align:center",
                                       h4("Submitting...")
                                  ),
                                  div(id = "error",
                                      style = "text-align:center",
                                      div(br(), tags$b("Error: "), span(id = "error_msg"))
                                  ),
                                  div(id = "wait_msg",
                                      style = "text-align:center",
                                      h4("Please, allow a few minutes to download your data."),
                                      hr() 
                                  ),
                                  div(id="goToNextTab",
                                      style = "text-align:center",
                                      h4(uiOutput( outputId = "url")), 
                                      h4("Your data has been downloaded. Please, go to the Data Description tab and fill in the form."),
                                      hr()
                                  ),
                                  div(id="falseGSE",
                                      style = "text-align:center",
                                      h3(textOutput(outputId = "falseGSE")),
                                      hr()
                                  ),
                                  div(id="data_sub_msg",
                                      style = "text-align:center",
                                      h4("Your data has been uploaded. Please, go to the Data Description tab and fill in the form."),
                                      hr()
                                  )
                                ),
                                bsTooltip("url","Click to redirect to GEO website",placement = "top", trigger = "hover")
                         ),
                         h4("You can access the source code and/or the docker image using the following links:", style="text-align: left"),
                         h4(uiOutput (outputId= "GitHub"), style="text-align:left"),
                         h4(uiOutput (outputId= "DockerHub"), style="text-align:left")
              ),
              
              tabPanel( title = "Data Description",
                        br(),
                        br(),
                        br(),
                        br(),
                        div(id= "firstInput1",
                            style = "text-align:center",
                            h4("You should input your data first."),
                            hr()
                        ),
                        shinyjs::hidden(
                          div(id="wait_msg4",
                              style = "text-align:center",
                              h4("Please, allow a few minutes to download your data."), 
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id= "choosing_platform",
                              uiOutput( outputId= "choosePlatform"),
                              tags$head(tags$style(HTML("#submit_platform{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;}"))),
                              actionButton(inputId = "submit_platform", label = "Submit"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id= "data_table",
                              dataTableOutput( outputId = "df"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id= "data_table_up",
                              DT::dataTableOutput( outputId = "df_up"),
                              helpText("Please, use valid names for treatment.", 
                                       "A syntactically valid name consists of letters, numbers, the dot or underscore characters,", 
                                       "and starts with a letter or the dot not followed by a number. Spaces are not allowed."),
                              helpText("Columns 'treatment' and 'replicate' should be filled out!"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="changes",
                              style = "text-align:center",
                              h4("Double click to complete the table and then press Save changes"),
                              tags$head(tags$style(HTML("#changes_done{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;}"))),
                              tags$head(tags$style(HTML("#discard{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;}"))),
                              actionButton( inputId = "changes_done", label = "Save changes"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="no_comparison",
                              style = "text-align:center",
                              h4(textOutput ( outputId = "comp")),
                              hr()
                          ),
                          shinyjs::hidden(
                            div(id="comparison",
                                uiOutput(outputId= "compBut"),
                                tags$head(tags$style(HTML("#submit_comp{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;}"))),
                                actionButton( inputId="submit_comp", label = "Submit"),
                                hr()
                            )
                          ),
                          shinyjs::hidden(
                            div(id="ctr",
                                uiOutput(outputId="comp_ctr"),
                                tags$head(tags$style(HTML("#submit_ctr{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;}"))),
                                actionButton( inputId="submit_ctr", label = "Submit"),
                                hr()
                            )
                          ),
                          shinyjs::hidden(
                            div(id="tr",
                                uiOutput(outputId="comp_tr"),
                                tags$head(tags$style(HTML("#submit_tr{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;}"))),
                                actionButton( inputId="submit_tr", label = "Submit"),
                                hr()
                            )
                          ),
                          shinyjs::hidden(
                            div(id="btns",
                                radioButtons( inputId= "adjMet", label = "Select the adjustment method", 
                                              choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                              selected = "fdr"),
                                div(id="slider_btns",
                                    tags$head(tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                                background: #B39DDB; border:0px;}
                                                                .irs-from, .irs-to, .irs-single { background: black }"))),
                                    sliderInput( inputId= "logFC", label = "Set the absolute logFC threshold:", min = 0.001, max = 2.0, step = 0.01,
                                                 value = 0.5),
                                    tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {
                                                        background: #B39DDB; border:0px;}
                                                      .irs-from, .irs-to, .irs-single { background: black }'")),
                                    sliderInput( inputId= "adjPval", label = "Set the adjusted P value:", min = 0.001, max = 0.1, step= 0.001, 
                                                 value = 0.05)
                                ),
                                hr()
                            ),
                            div(id= "run",
                                style = "text-align:center",
                                h4("press Run the analysis"),
                                tags$head(tags$style(HTML('#anls{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;}'))),
                                actionButton( inputId= "anls", label = "Run the analysis"),
                                hr()
                            ),
                            bsPopover("adjMet",
                                      title="<b>Adjustment for multiple comparisons</b>",
                                      content=c("The probability of committing false statistical inferences",
                                                "would considerably increase when more than one hypothesis",
                                                "are simultaneously tested (<b>multiple comparisons</b>), in which", 
                                                "case proper adjustment is required. Adjustment methods include", 
                                                "the Bonferroni correction (<code>bonferroni</code>) in which p-values are", 
                                                "multiplied by the number of comparisons. Less conservative", 
                                                "corrections are also included by Holm (1979) (<code>holm</code>),", 
                                                "Hochberg (1988) (<code>hochberg</code>), Hommel (1988) (<code>hommel</code>), ",
                                                "Benjamini & Hochberg (1995) (<code>BH</code> or its alias <code>fdr</code>), and ", 
                                                "Benjamini & Yekutieli (2001) (<code>BY</code>). A pass-through option ",
                                                "(<code>none</code>) is also included. For a more detailed review, see: ",
                                                "<b>Shaffer, J. P. (1995). Multiple hypothesis testing. ",
                                                "<i>Annual Review of Psychology,</i> 46, 561--576.</b>"),
                                      placement="right",
                                      trigger = "hover"),
                            bsPopover("logFC",
                                      title = "<b>Absolute log2 Fold Change</b>",
                                      content = c("Fold change is used to measure change on the expression ",
                                                  "level of a gene. Here, we use the fold change on logarithmic scale ",
                                                  "(base=2). You can set the absolute value of logFC so that the gene ",
                                                  "expression is considered differential as compared to the control."),
                                      placement = "right", 
                                      trigger = "hover"),
                            bsPopover("adjPval",
                                      title = "<b>Adjusted P Value</b>",
                                      content = c("You can set the adjusted p-value threshold in order ",
                                                  "to define statistical significance."),
                                      placement = "right", 
                                      trigger = "hover"),
                            
                            shinyjs::hidden(
                              div(id="wait_msg2",
                                  style = "text-align:center",
                                  h4("Running..."),
                                  h4("Please, wait a few minutes."),
                                  hr()
                              )
                            ),
                            shinyjs::hidden(
                              div(id="msg",
                                  style = "text-align:center",
                                  h4("The analysis is done. Please, go to the Results tab."),
                                  hr()
                                  )
                              )
                            )
                          )
                        ),
              tabPanel(title = "Results",
                                    br(),
                                    br(),
                                    br(),
                                    br(),
                                    div(id= "firstInput2",
                                        style = "text-align:center",
                                        h4("You should input your data first."),
                                        hr()
                                        ),
                                    shinyjs::hidden(
                                      div(id= "firstRun",
                                          style = "text-align:center",
                                          h4("Complete the Data Description form and then press Run the analysis."),
                                          hr()
                                          )
                                      ),
                                    shinyjs::hidden(
                                      div(id="wait_msg3",
                                          style = "text-align:center",
                                          h4("Please, wait a few minutes until the analysis is done."), 
                                          hr()
                                          )
                                      ),
                                    shinyjs::hidden(
                                      div(id="noDEGs",
                                          style = "text-align:center",
                                          h3("There are NO differentially expressed genes."), 
                                          hr()
                                          )
                                      ),
                                    shinyjs::hidden(
                                      div(id="no_annotation",
                                          style = "text-align:center",
                                          h3("Annotation file not found"), 
                                          hr()
                                          )
                                      ),
                       shinyjs::hidden(
                         div(
                           id = "download",
                           fluidRow(
                                    style = "text-align:center;",
                                    h3("The results can be downloaded as:",
                                       HTML("<b>visualization plots</b>"),
                                       " (including histogram of adjusted p value against the number of probes, boxplot, 
                                       interactive heatmap, interactive volcano plot, and PCA plots, i.e., 
                                       the scree plot, the grouping of samples against PC1 and PC2, and the biplot),", 
                                       HTML("the DEGs list as a "),
                                       HTML("<b>.csv file,</b>"),
                                       HTML("or the DEGs list as a "),
                                       HTML("<b>.tsv file</b>"),
                                       HTML("by pressing the corresponding button.")
                                       ),
                                    br(), br(),
                                    tags$head(tags$style(HTML('#down_plots{background-color:#B39DDB; color:#000000;font-weight:bold; border:0px;}'))),
                                    downloadButton(outputId = "down_plots", label = "Visualization plots"),
                                    tags$head(tags$style(HTML('#down_csv{background-color:#B39DDB; color:#000000;font-weight:bold; border:0px;}'))),
                                    downloadButton(outputId = "down_csv", label = ".csv file"),
                                    tags$head(tags$style(HTML('#down_tsv{background-color:#B39DDB; color:#000000;font-weight:bold; border:0px;}'))),
                                    downloadButton(outputId = "down_tsv", label = ".tsv file"),
                                    tags$head(tags$style(HTML('#table_info{background-color:#B39DDB; color:#000000;font-weight:bold; border:0px;}'))),
                                    bsButton(inputId = "table_info", 
                                             label = icon(name = "info", class = ".fa-info"),
                                             style = "info",
                                             type = "toggle",
                                             block = FALSE
                                             ), 
                                    hr()
                                    )
                             )
                           ),
                                    shinyjs::hidden(
                                      div(id="DEGs",
                                          DT::dataTableOutput( outputId = "DEGs"),
                                          hr()
                                          )
                                      ),
                       bsTooltip("down_csv", "Click to download the DE genes list as a comma separated file", 
                                 placement = "left", trigger = "hover"),
                       bsTooltip("down_tsv", "Click to download the DE genes list as a tab separated file", 
                                 placement = "left", trigger = "hover"),
                       bsTooltip("down_plots", "Click to download the visualization plots as a .zip file", 
                                 placement = "left", trigger = "hover"),
                       bsPopover (id= "table_info",
                                  title = "<b>Table Info</b>",
                                  content =c("<p><b>probeID</b>= Affymetrix unique identifier assigned ",
                                             "to each probe on the array</p> <p><b>gene symbol</b>= the official ",
                                             "gene symbol name used in Genomics databases</p> <p><b>logFC</b>= ",
                                             "the value of the contrast; this represents a log2-fold ",
                                             "change between two experimental conditions, i.e. control ",
                                             "and treated sample</p> <p><b>AveExpr</b>= the average ",
                                             "log2-expression level for that gene across all the arrays ",
                                             "in the experiment</p> <p><b>t</b>= the moderated t-statistic</p> ",
                                             "<p><b>P.Value</b>= the associated p-value</p>  <p><b>adj.P.Value</b>= ",
                                             "the p-value adjusted for multiple testing</p>",
                                             "<p><b>B</b>= log-odds that the gene is differentially expressed</p>"),
                                  placement = "bottom", 
                                  trigger = "click")
                       ),
              tabPanel( title= "WebGestalt Over-Representation Analysis",
                        br(),
                        br(),
                        br(),
                        br(),
                        shinyjs::hidden(
                          div(id="wait_msg5",
                              style = "text-align:center",
                              h4("You should find the differentially expressed genes before you proceed to WebGestalt ORA"), 
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="WebGestaltR_but",
                              style = "text-align:center",
                              tags$head(tags$style(HTML('#WebGestaltR{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;width:200px;margin:auto;text-align:center;}'))),
                              actionButton( inputId= "WebGestaltR", label = "WebGestalt ORA"),
                              bsTooltip(id= "WebGestaltR", 
                                        title = "Press to proceed to WebGestalt Over-Representation Analysis (ORA)",
                                        placement = "bottom",
                                        trigger = "hover"),
                              br(),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id= "choosing_organism", 
                              uiOutput( outputId= "chooseOrganism"),
                              tags$head(tags$style(HTML("#submit_Organism{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;}"))),
                              actionButton(inputId = "submit_Organism", label = "Submit"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id= "choosing_RefSet",
                              uiOutput( outputId= "chooseRefSet"),
                              tags$head(tags$style(HTML("#submit_RefSet{background-color:#B39DDB; color:#000000;font-weight:bold;border:0px;margin:auto;}"))),
                              actionButton(inputId = "submit_RefSet", label = "Submit"),
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="wait_msg6",
                              style = "text-align:center",
                              h4("Please, wait a few minutes until the WebGestalt analysis is done."), 
                              hr()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="WebGestalt_down",
                              style = "text-align:center",
                              tags$head(tags$style(HTML('#down_WebGestalt{background-color:#B39DDB; color:#000000;font-weight:bold; border:0px;}'))),
                              downloadButton( outputId = "down_WebGestalt", label = "Download WebGestalt ORA results"),
                              br(),
                              br(),
                              br()
                          )
                        ),
                        shinyjs::hidden(
                          div(id="WebGestalt_res",
                              uiOutput( outputId = "WebGestalt")
                          )
                        ),
                        shinyjs::hidden(
                          div(id="msg_no_annot",
                              style = "text-align:center",
                              h4("You cannot proceed to WebGestalt ORA since there is no annotation for your data"), 
                              hr()
                          )
                        )
              ),
              tabPanel( title = "About",
                        br(),
                        br(),
                        br(),
                        br(),
                        div ( id= "about",
                              uiOutput(outputId= "about")
                        ),
                        hr(),
                        div ( id="email",
                              style = "text-align:left",
                              h4("Email us: ",a ("dexplore.app[at]gmail.com", 
                                                 href="mailto:dexplore.app[at]gmail.com", target="_blank"))
                        ),
                        hr(),
                        div ( id = "users_guide",
                              style = "text-align:center",
                              uiOutput(outputId= "users_guide_msg"),
                              uiOutput(outputId= "users_guide_view")
                              #tags$head(tags$style(HTML('#users_guide{background-color:#B39DDB; color:#000000; 
                              #                           font-weight:bold; border:0px; width:125px; padding:3px 0; text-align:center; margin:auto;}'))),
                              #downloadButton( outputId="users_guide", label = "User's guide")
                              #tags$iframe(style="height:400px; width:100%; scrolling=yes",
                              #             src="c:/Users/Anna/Desktop/PhD/git/DExplore/GSE/Users_guide.pdf")
                        ),
                        hr()
                        )
              )
  )
)
