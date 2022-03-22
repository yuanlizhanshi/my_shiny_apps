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
library(ggprism)
library(ggpubr)
determine_p <- function(pvalue){
  case_when(
    pvalue <= 0.05 & pvalue > 0.01 ~ paste0("pvalue =",pvalue,'   ',"*"),
    pvalue <= 0.01 & pvalue > 0.001~ paste0("pvalue =",pvalue,'   ',"**"),
    pvalue <= 0.001 ~ paste0("pvalue =",pvalue,'  ',"***"),
    TRUE ~ paste0("pvalue =",pvalue,'   ',"NS")
  ) 
}
bar_plot <- function(df,pvalue,x_labels = c('Control','Treat')){
  pvalues <- determine_p(pvalue)
  ggplot(df,aes(sample,Fold,fill= sample)) + 
    geom_col(width = 0.5,position = position_dodge(width = 0.1)) + theme_classic() +
    geom_errorbar(aes(ymin = Fold - sd, ymax = Fold + sd),width=0.2) +
    labs(x = '',y = 'Relative expression levels')+
    guides(fill = 'none')+
    theme_prism()+
    scale_fill_manual(values = c("#69b3a2","#836FFF"))+
    scale_x_discrete(labels = x_labels)+
    ylim(NA,max(df$Fold,na.rm = T)+0.2)+
    theme(axis.title.x = element_blank(),
          aspect.ratio = 1.3)+
    geom_signif(comparisons = list(c("Control_group", "Treat_group")),
                annotations = pvalues,
                margin_top = 0.3,
                tip_length = 0.1,
                vjust = -0.7
    )
}

ui <- fluidPage(
  div(
    fluidRow(column(2,textInput("gapdh", "内参基因(输入内参基因的名字)",'GAPHD',width = '200px'))),
    fluidRow(
      column(1,numericInput("ref_gapdh1", "对照组内参基因Cq1:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gapdh2", "对照组内参基因Cq2:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gapdh3", "对照组内参基因Cq3:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gapdh4", "对照组内参基因Cq4:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gene1", "对照组目的基因Cq1:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gene2", "对照组目的基因Cq2:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gene3", "对照组目的基因Cq3:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("ref_gene4", "对照组目的基因Cq4:", NA, min = 1, max = 100,width = '80px'))
    ),
    fluidRow(column(2,textInput("targetgene", "目的基因(输入目的基因的名字)",NA,width = '200px'))),
    fluidRow(
      column(1,numericInput("treat_gapdh1", "处理组内参基因Cq1:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gapdh2", "处理组内参基因Cq2:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gapdh3", "处理组内参基因Cq3:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gapdh4", "处理组内参基因Cq4:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gene1", "处理组目的基因Cq1:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gene2", "处理组目的基因Cq2:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gene3", "处理组目的基因Cq3:", NA, min = 1, max = 100,width = '80px')),
      column(1,numericInput("treat_gene4", "处理组目的基因Cq4:", NA, min = 1, max = 100,width = '80px'))
    ),
    style="font-size:80%"
  ),
  tableOutput("values"),
  "结果",
  verbatimTextOutput("pvalue"),
  tableOutput("results_df1"),
  fluidRow(
    helpText('滑动调节图片长宽'),
    column(6,sliderInput("height", "height", min = 100, max = 1000, value = 300)),
    column(6,sliderInput("width", "width", min = 100, max = 1000, value = 300))
  ),
  fluidRow(
    column(2,textInput("control_name", "输入对照样本的名字","Control",width = '200px')),
    column(2,textInput("treat_name", "输入处理样本的名字","Treat",width = '200px'))
  ),
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

server <- function(input, output) {
  target_id <- reactive({input$targetgene})
  rowcq_df_obj <- reactive({
    gadph_cqvalues <- na.omit(c(input$ref_gapdh1,input$ref_gapdh2,input$ref_gapdh3,input$ref_gapdh4,
                                input$treat_gapdh1,input$treat_gapdh2,input$treat_gapdh3,input$treat_gapdh4))
    target_cqvalues <- na.omit(c(input$ref_gene1,input$ref_gene2,input$ref_gene3,input$ref_gene4,
                                 input$treat_gene1,input$treat_gene2,input$treat_gene3,input$treat_gene4))
    ref_id <- input$gapdh
    target_id <- input$targetgene
    reps <- length(c(gadph_cqvalues,target_cqvalues))/4
    samples <- c(rep('Control_group',reps),rep('Treat_group',reps))
    rowcq_df <- data.frame(
      sample = samples,
      cqGAPHD = gadph_cqvalues,
      cqgene = target_cqvalues
    )
    colnames(rowcq_df)[2] <- paste0(ref_id,'_Cq')
    colnames(rowcq_df)[3] <- paste0(target_id,'_Cq')
    rowcq_df$mean_cq1 <- rowcq_df[[3]]
    rowcq_df$mean_cq2 <- rowcq_df[[2]]
    rowcq_df <- rowcq_df %>%
      mutate(deltaCT = mean_cq1-mean_cq2) %>%
      mutate(two_NegdeltaCT = 2**(-deltaCT)) %>%
      mutate(mean_two_NegdeltaCT = c(rep(mean(two_NegdeltaCT[1:reps]),reps),
                                     rep(mean(two_NegdeltaCT[c(reps+1):c(reps*2)]),reps))) %>%
      mutate(two_NegdeltadeltaCT = two_NegdeltaCT/mean_two_NegdeltaCT[1]) 
    return(rowcq_df)
  })
  fold_res <- reactive({
    fold_res_df <- rowcq_df_obj() %>% group_by(sample) %>% 
      summarise(Fold = mean(two_NegdeltadeltaCT),sd = sd(two_NegdeltadeltaCT))
    if(fold_res_df$Fold[2] <= 1){
      fold_res_df$Fold <- fold_res_df$Fold/fold_res_df$Fold[2]
    }
    return(fold_res_df)
  })
  pvalue <- reactive({
    reps <- nrow(rowcq_df_obj())/2
    ttest_res <- t.test(
      rowcq_df_obj()[1:reps,ncol(rowcq_df_obj())],
      rowcq_df_obj()[c(reps+1):c(reps*2),ncol(rowcq_df_obj())]
    )
    return(round(ttest_res$p.value,4))
  })
  x_label <- reactive({
    return(c(input$control_name,input$treat_name))
  })
  output$values <- renderTable(
    {
      rowcq_df_obj()
    }
    ,digits = 4)
  output$results_df1 <- renderTable(
    {
      fold_res()  
    }
    ,digits = 4)
  output$pvalue <- renderText({
    determine_p(pvalue())
  }
  )
  output$barplots <- renderPlot(
    width = function() input$width,
    height = function() input$height,
    {
      plot(bar_plot(fold_res(),pvalue(),x_labels = x_label()))
    },res = 96)
  output$downloadPlot <- downloadHandler(
    filename = function(){paste0(target_id(),'.tiff')},
    content = function(file) {
      tiff(file,width = 3000, height = 4000, units = "px", res = 600, compression = "lzw")
      plot(bar_plot(fold_res(),pvalue()))
      dev.off()
    }) 
}
shinyApp(ui, server)
