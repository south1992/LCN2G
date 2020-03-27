library(shiny)
library(shinythemes)
library(shinyjs)
library(visNetwork)
library(DT)
source("utils.R")

navbarPage("LC-N2G",
           shinyjs::useShinyjs(),
           theme = shinytheme("cerulean"),
           tabPanel("Data Preprocess",
                    sidebarLayout(
                      sidebarPanel(
                        numericInput("mcpm","Filter out gene expression(cpm) below:",5),
                        numericInput("maxcpm","Filter out gene expression(cpm) above:",500),
                        numericInput("mvar","Filter out gene expression(cpm) sd below :",0.1),
                        radioButtons("mnorm","Normalization method",choices = c("cpm" = "cpm","log2 cpm" = "lcpm"),selected = "cpm"),
                        sliderInput("mabscor","Filter out gene expression(cpm) maxium absolute correlation with nutrition variables below:",min = 0,max = 1,step = 0.05,value = 0.6),
                        actionButton("preprocess", "Analysis")
                      ),
                      mainPanel(
                        plotOutput("Preprocess",height = "700px"),
                        textOutput("ngene")
                      )
                    )
                    ),
           
           tabPanel("Gene Cluster",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText("Paramerters for clustering using gene-gene relationship:"),
                          radioButtons("mclust1","Cluster method:",
                                       c("Cluster by gene expression" = "m0",
                                         "Cluster by gene-gene correlation(WGCAN) " = "m1"),selected = "m1")
                        ),
                        wellPanel(
                          helpText("Parameters for clustering using gene-nutrition relationship:"),
                          radioButtons("mclust2","Cluster method:",
                                       c("Cluster by gene-nutrition covariance " = "m0",
                                         "Cluster by gene-nutrition correlation" = "m1"),selected = "m1")
                        ),
                        wellPanel(
                          helpText("Other Parameters for clustering"),
                          sliderInput("alpha","mix coefficient for two kind of cluster method(0 for total gene-gene,1 for total gene-nutrition)",
                                      min = 0,max = 1,step = 0.05,value = 0.5),
                          radioButtons("scaledist","Scaling method for distance matrix:",
                                       c("No scaling" = "m0",
                                         "Linear Transformation scaling" = "m1",
                                         "Exponential Transformation scaling" = "m2"),selected = "m2"),
                          radioButtons("mcut","Tree Cutting Method:",
                                       c("Dynamic tree cut" = "c1",
                                         "By threshold" = "c2")),
                          textInput("k","Desired number of groups",value = 10),
                          sliderInput("thresh","threshold of correlation for network visualization for clustering result:",
                                      min = 0,max = 1,step = 0.05,value = 0.5),
                          actionButton("cluster", "Analysis")
                        )
                      ),
                      mainPanel(
                        tabsetPanel(
                          id = "Clust_output",
                          tabPanel("Dendro",plotOutput("denro"),DT::dataTableOutput("clusttab")),
                          tabPanel("Table",DT::dataTableOutput("exprtab")),
                          tabPanel("Summary",visNetworkOutput("vnet",height = "700px"))
                        )
                        )
                    )),
           
           tabPanel("Gene Nutrition Visualization",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                          helpText("Parameters for zaxis"),
                          radioButtons("zinput","z-axis choose by:",
                                       c("Cluster" = "z1","Gene name" = "z2")),
                          uiOutput("zaxis10"),
                          textInput("zaxis2","z-axis(gene):",value = "Ucp2"),
                          checkboxInput("scalez","zscore for z-axis:",value = T)
                        ),
                        
                        wellPanel(
                          helpText("Parameters for minimum local consistency nutrition variables selection"),
                          selectInput("Nutrivar","Nutrition variables considered",
                                              choices = choices_xy(),selected = choices_xy0(),multiple = T),
                          uiOutput("m0"),
                          actionButton("visualization0", "Analysis")
                                  ),
                        
                        wellPanel(
                          helpText("Parameters for 2D visualization of Nutrigenetics"),
                          uiOutput("xaxis0"),
                          uiOutput("yaxis0"),
                          checkboxInput("scalexy","zscore for xy-axis:",value = T),
                          actionButton("visualization", "Analysis")
                        )
                        
                      ),
                      mainPanel(
                        textOutput("opt"),
                        plotOutput("nutrivisual",height = "700px"))
                    ))
  
)