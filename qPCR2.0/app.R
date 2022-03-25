#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggplot2)

get_list <- function(gene_name,df = all_gene_cq ){
  df <- data.frame(
    sample = df[[1]],
    GAPDH = df[[2]],
    target_gene = df[[gene_name]]
  )
  return(df)
}

single_calcluate2 <- function(dfx){
  repeat_num <- nrow(dfx)/2
  df1 <- dfx %>% pivot_longer(cols = -sample,names_to = 'gene',values_to = 'cq') %>% group_by(gene) %>% mutate(n = row_number())
  df2 <- df1 %>% pivot_wider(names_from = c(gene),values_from = cq,names_prefix = 'cq')
  df3 <- df2 %>% mutate(mean_cq1 = unlist(c(df2[str_which(colnames(df2),'target_gene')],use.names=F)),
                        mean_cq2 = unlist(c(df2[str_which(colnames(df2),'GAPDH')],use.names=F))) %>%
    mutate(deltaCT = mean_cq1 - mean_cq2) %>%
    mutate(two_NegdeltaCT = 2**(-deltaCT)) %>%
    mutate(mean_two_NegdeltaCT = c(rep(mean(two_NegdeltaCT[1:repeat_num]),repeat_num),
                                   rep(mean(two_NegdeltaCT[c(repeat_num+1):c(repeat_num*2)]),repeat_num))) %>%
    mutate(two_NegdeltadeltaCT = two_NegdeltaCT/mean_two_NegdeltaCT[1])
  return(df3)
}
plot_df <- function(df3){
  determine_p <- function(pvalue){
    case_when(
      pvalue <= 0.05 & pvalue > 0.01 ~ "*",
      pvalue <= 0.01 & pvalue > 0.001~ "**",
      pvalue <= 0.001 ~ "***",
      TRUE ~ "NS"
    ) 
  }
  repeat_num = nrow(df3)/2
  names <- unique(df3[[1]])
  sample_num = length(names)
  ttest_res <- t.test(df3[1:repeat_num,10],df3[c(repeat_num+1) : c(repeat_num*2),10])
  aa <- data.frame(
    sample = names,
    Fold = c(
      apply(df3[1:repeat_num,10],2, mean),
      apply(df3[c(repeat_num+1) : c(repeat_num*2),10],2, mean)
    ),
    SD =   c(
      apply(df3[1:repeat_num,10],2, sd),
      apply(df3[c(repeat_num +1):c(repeat_num*2),10],2, sd)),
    pvalue = round(ttest_res$p.value,4),
    significant = determine_p(round(ttest_res$p.value,4))
  )
  return(aa)
}
bar_plot <- function(test_res){
  get_maxfold <- function(test_res){
    df1 <- test_res %>% group_by(gene_name) %>% summarise(max_fold = max(Fold))
    df2 <- test_res %>% group_by(gene_name) %>% summarise(max_SD = max(SD))
    return(df1$max_fold+df2$max_SD +0.2)
  }
  get_ano <- function(test_res){
    df1 <- test_res %>% group_by(gene_name) %>% summarise(ano = unique(significant))
    return(df1$ano)
  }
  ggplot(test_res,aes(gene_name,Fold,fill= sample)) + 
    geom_bar(stat = "identity",position = position_dodge(width = 1)) + theme_classic() +
    geom_errorbar(aes(ymin = Fold - SD, ymax = Fold + SD),
                  position=position_dodge(width=1), 
                  width=0.3,size=0.3,colour="black") +
    labs(x = '',y = 'Relative expression levels')+
    guides()+
    ggprism::theme_prism()+
    scale_fill_manual(values = c("#69b3a2","#836FFF"))+
    theme(axis.title.x = element_blank(),
          aspect.ratio = 1.3,
          axis.text.x = element_text(angle = 45))+
    ggpubr::geom_signif(y_position = get_maxfold(test_res), 
                        xmin = rep(1:length(get_ano(test_res))) -0.2, 
                        xmax = rep(1:length(get_ano(test_res))) + 0.2,
                        annotation = get_ano(test_res),
                        tip_length = 0)
}

ui <- fluidPage(
  
  # Application title
  titlePanel("qRT-PCR全自动画图"),
  fluidRow(
      column(6,
             fileInput("upload", NULL, buttonLabel = "Upload...", multiple = FALSE,accept = c(".xls", ".xlsx"))
             ),
      column(3,
             downloadButton("download_test", "Download excel format")
             )
    ),
    div(
      fluidRow(
        column(6,
               dataTableOutput("preview")
        )
      ),
      style = "font-size:80%"
  ),
  div(
    fluidRow(
      "Results"
      ),
    style = "font-size:200%"
  ),
  tableOutput("results"),
  fluidRow(
    plotOutput("barplots")
  ),
  fluidRow(
    downloadButton('downloadPlot', 'Download Plot')
  ),
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           xls = readxl::read_xls(input$upload$datapath),
           xlsx = readxl::read_xlsx(input$upload$datapath),
           validate("Invalid file; Please upload a .xls or .xlsx file")
    )
  })
  all_gene_list <- reactive({
    req(input$upload)
    map(colnames(data())[3:ncol(data())],get_list,df = data()) %>% setNames(colnames(data())[3:ncol(data())])
  })
  calculation_process <- reactive({
      req(input$upload)
      map_dfr(all_gene_list(),single_calcluate2,.id = 'gene_name')
    })
  calculation_results <- reactive({
    req(input$upload)
    map_dfr(map(all_gene_list(),single_calcluate2),plot_df,.id = 'gene_name')
  })
  output$preview <- renderDataTable({
    calculation_process()
  },
  options = list(pageLength = 10)
  )
  output$results <- renderTable(
    {
    calculation_results()
    }
    ,digits = 4)
  output$barplots <- renderPlot({
    bar_plot(calculation_results())
  })
  output$downloadPlot <- downloadHandler(
    filename = function(){paste0(paste0(colnames(data())[c(-1,-2)],collapse = '_'),'.tiff')},
    content = function(file) {
      tiff(file,width = 3000, height = 4000, units = "px", res = 600, compression = "lzw")
      plot(bar_plot(calculation_results()))
      dev.off()
    }) 
  output$download_test <- downloadHandler(
    filename = function(){
      "format.xlsx"
    },
    content = function(file){
      file.copy("format.xlsx", file)
    }
  )
}


# Run the application 
shinyApp(ui = ui, server = server)