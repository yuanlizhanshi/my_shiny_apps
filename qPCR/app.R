library(shiny)
library(tidyverse)
library(shinydashboard)

get_list <- function(gene_name,df = all_gene_cq ){
  df <- data.frame(
    sample = df[[1]],
    GAPDH = df[[2]],
    target_gene = df[[gene_name]]
  )
  return(df)
}
single_calcluate2 <- function(dfx){
  repeat_num <- nrow(dfx)/length(table(dfx[[1]]))
  df1 <- dfx %>% pivot_longer(cols = -sample,names_to = 'gene',values_to = 'cq') %>% group_by(gene) %>% mutate(n = row_number())
  df2 <- df1 %>% pivot_wider(names_from = c(gene),values_from = cq,names_prefix = 'cq')
  df3 <- df2 %>% mutate(mean_cq1 = unlist(c(df2[str_which(colnames(df2),'target_gene')],use.names=F)),
                        mean_cq2 = unlist(c(df2[str_which(colnames(df2),'GAPDH')],use.names=F))) %>%
    mutate(deltaCT = mean_cq1 - mean_cq2) %>%
    mutate(two_NegdeltaCT = 2**(-deltaCT))
  sample_mean_two_NegdeltaCT <- df3 %>% group_by(sample) %>% summarise(mean_two_NegdeltaCT = mean(two_NegdeltaCT))
  df5 <- df3 %>% mutate(mean_two_NegdeltaCT = unlist(map(sample_mean_two_NegdeltaCT$mean_two_NegdeltaCT,rep,repeat_num)))  %>%
    mutate(two_NegdeltadeltaCT = two_NegdeltaCT/mean_two_NegdeltaCT[1])
  return(df5)
}
plot_df_1 <- function(df3){
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
bar_plot_1 <- function(test_res){
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



plot_df_2 <- function(df3){
  fold_value <- df3 %>% group_by(sample) %>% summarise(Fold = mean(two_NegdeltadeltaCT),SD = sd(two_NegdeltadeltaCT))
  return(fold_value)
}
bar_plot_2 <- function(test_res){
  ggplot(test_res,aes(gene_name,Fold,fill= sample)) + 
    geom_bar(stat = "identity",position = position_dodge(width = 1)) + theme_classic() +
    geom_errorbar(aes(ymin = Fold - SD, ymax = Fold + SD),
                  position=position_dodge(width=1), 
                  width=0.3,size=0.3,colour="black") +
    labs(x = '',y = 'Relative expression levels')+
    theme(axis.title.x = element_blank(),
          aspect.ratio = 1.3,
          axis.text.x = element_text(angle = 45))+
    ggprism::theme_prism()
}



header <- dashboardHeader(
  title = "qPCR Calculator"
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Paired sample", tabName = "Paired_sample"),
    menuItem("Multi sample", tabName = "Multi_sample"),
    div(style="text-align:center")
  )
  
)

body <- dashboardBody(
  ##change the CSS format of sidebar 
  tags$head( 
    tags$style(HTML(".main-sidebar { font-size: 20px };")) #change the font size to 20
  ),
  
  tabItems(
    # Paired_sample tab content
    
    tabItem(tabName = "Paired_sample",
            fluidRow(
              column(6,
                     fileInput("upload_1", NULL, buttonLabel = "Upload Paired sample...", multiple = FALSE,accept = c(".xls", ".xlsx"))
              ),
              column(3,
                     downloadButton("download_test_1", "Download excel format")
              )
            ),
            div(
              fluidRow(
                column(6,
                       dataTableOutput("preview_1")
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
            tableOutput("results_1"),
            fluidRow(
              plotOutput("barplots_1")
            ),
            fluidRow(
              downloadButton('downloadPlot', 'Download Plot')
            ),
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
            )
    ),
    
    
    # Multi_sample tab content
    tabItem(tabName = "Multi_sample",
            fluidRow(
              column(6,
                     fileInput("upload_2", NULL, buttonLabel = "Upload Multi sample...", multiple = FALSE,accept = c(".xls", ".xlsx"))
              ),
              column(3,
                     downloadButton("download_test_2", "Download excel format")
              )
            ),
            div(
              fluidRow(
                column(6,
                       dataTableOutput("preview_2")
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
            tableOutput("results_2"),
            div(
              fluidRow(
                "Pairwise t test"
              ),
              style = "font-size:200%"
            ),
            tableOutput("ttest"),
            fluidRow(
              plotOutput("barplots_2")
            ),
            fluidRow(
              downloadButton('downloadPlot_2', 'Download Plot')
            ),
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
            )
    )
  )
)

ui <-  dashboardPage(header, sidebar, body,
                     skin = "purple")

server <- function(input, output) {
  data <- reactive({
    req(input$upload_1)
    ext <- tools::file_ext(input$upload_1$name)
    switch(ext,
           xls = readxl::read_xls(input$upload_1$datapath),
           xlsx = readxl::read_xlsx(input$upload_1$datapath),
           validate("Invalid file; Please upload a .xls or .xlsx file")
    )
  })
  all_gene_list <- reactive({
    req(input$upload_1)
    map(colnames(data())[3:ncol(data())],get_list,df = data()) %>% setNames(colnames(data())[3:ncol(data())])
  })
  calculation_process <- reactive({
    req(input$upload_1)
    map_dfr(all_gene_list(),single_calcluate2,.id = 'gene_name')
  })
  calculation_results <- reactive({
    req(input$upload_1)
    map_dfr(map(all_gene_list(),single_calcluate2),plot_df_1,.id = 'gene_name')
  })
  output$preview_1 <- renderDataTable({
    calculation_process()
  },
  options = list(pageLength = 10)
  )
  output$results_1 <- renderTable(
    {
      calculation_results()
    }
    ,digits = 4)
  output$barplots_1 <- renderPlot({
    bar_plot_1(calculation_results())
  })
  
  ####server funcuoint of paired 
  data2 <- reactive({
    req(input$upload_2)
    ext <- tools::file_ext(input$upload_2$name)
    switch(ext,
           xls = readxl::read_xls(input$upload_2$datapath),
           xlsx = readxl::read_xlsx(input$upload_2$datapath),
           validate("Invalid file; Please upload a .xls or .xlsx file")
    )
  })
  all_gene_list2 <- reactive({
    req(input$upload_2)
    map(colnames(data2())[3:ncol(data2())],get_list,df = data2()) %>% 
      setNames(colnames(data2())[3:ncol(data2())])
  })
  calculation_process2 <- reactive({
    req(input$upload_2)
    map_dfr(all_gene_list2(),single_calcluate2,.id = 'gene_name')
  })
  calculation_results2 <- reactive({
    req(input$upload_2)
    map_dfr(map(all_gene_list2(),single_calcluate2),plot_df_2,.id = 'gene_name')
  })
  output$preview <- renderDataTable({
    calculation_results2()
  },
  options = list(pageLength = 10)
  )
  p_value <- reactive({
    calculation_process2() %>% rstatix::pairwise_t_test(two_NegdeltadeltaCT ~ sample, p.adjust.method = "bonferroni") %>% select(2,3,6,7)
  })
  output$preview_2 <- renderDataTable({
    calculation_process2()
  },
  options = list(pageLength = 10)
  )
  output$results_2 <- renderTable(
    {
      calculation_results2()
    }
    ,digits = 4)
  output$ttest <-renderTable(
    p_value()
  )
  output$barplots_2 <- renderPlot({
    bar_plot_2(calculation_results2())
  })
  
  ###download some thing
  output$downloadPlot <- downloadHandler(
    filename = function(){paste0(paste0(colnames(data())[c(-1,-2)],collapse = '_'),'.tiff')},
    content = function(file) {
      tiff(file,width = 3000, height = 4000, units = "px", res = 600, compression = "lzw")
      plot(bar_plot_1(calculation_results()))
      dev.off()
    })
  output$downloadPlot_2 <- downloadHandler(
    filename = function(){paste0(paste0(colnames(data())[c(-1,-2)],collapse = '_'),'.tiff')},
    content = function(file) {
      tiff(file,width = 3000, height = 4000, units = "px", res = 600, compression = "lzw")
      plot(bar_plot_2(calculation_results2()))
      dev.off()
    }) 
  output$download_test_1 <- downloadHandler(
    filename = function(){
      "format_1.xlsx"
    },
    content = function(file){
      file.copy("format_1.xlsx", file)
    }
  )
  output$download_test_2 <- downloadHandler(
    filename = function(){
      "format_2.xlsx"
    },
    content = function(file){
      file.copy("format_2.xlsx", file)
    }
  )
}

shinyApp(ui, server)
