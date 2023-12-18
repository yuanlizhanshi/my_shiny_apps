library(shiny)
library(tidyverse)
library(shinydashboard)
library(shinyWidgets)
library(dashboardthemes)
library(DT)
library(ggpubr)
plot_rna_exp <- function(gene){
    df <- NB4_gene_exp %>% filter(Gene == gene) %>% 
        pivot_longer(cols = -Gene,names_to = 'Stage',values_to = 'TPM') %>%
        mutate(day = str_extract(Stage,'A\\d')) %>%
        group_by(day) %>%
        mutate(mean_TPM = mean(TPM),
               sem = plotrix::std.error(TPM))
    ggplot(df) +
        geom_point(aes(day,TPM))+
        labs(x = NULL) + 
        geom_errorbar(aes(day,ymin = mean_TPM - sem, ymax = mean_TPM + sem), width=0.1) +
        geom_line(aes(day,mean_TPM,group=1)) +
        theme_bw()+
        theme(axis.text.x = element_text(size =15),
              axis.text.y = element_text(size =10),
              axis.title.y = element_text(size =15)) 
}
plot_prot_exp <- function(gene){
    df <- NB4_prot_exp %>% filter(Protein == gene) %>% 
        pivot_longer(cols = -Protein,names_to = 'Stage',values_to = 'Value') %>%
        mutate(day = str_extract(Stage,'NB4_A\\d')) %>%
        group_by(day) %>%
        mutate(mean_TPM = mean(Value),
               sem = plotrix::std.error(Value))
    
    ggplot(df) +
        geom_point(aes(day,Value))+
        labs(x = NULL) + 
        geom_errorbar(aes(day,ymin = mean_TPM - sem, ymax = mean_TPM + sem), width=0.1) +
        geom_line(aes(day,mean_TPM,group=1)) +
        theme_bw()+
        theme(axis.text.x = element_text(size =15),
              axis.text.y = element_text(size =10),
              axis.title.y = element_text(size =15)) 
}
plot_cor <- function(gene,method = 'pearson|spearman'){
    rna <- NB4_gene_exp %>% filter(Gene == gene) %>% 
        pivot_longer(cols = -Gene,names_to = 'Stage',values_to = 'TPM') %>%
        mutate(day = str_extract(Stage,'A\\d')) %>%
        group_by(day) %>%
        summarise(
            mean_TPM = mean(TPM))
    
    prot <- NB4_prot_exp %>% filter(Protein == gene) %>% 
        pivot_longer(cols = -Protein,names_to = 'Stage',values_to = 'Value') %>%
        mutate(day = str_extract(Stage,'NB4_A\\d')) %>%
        group_by(day) %>%
        summarise(
            mean_TPM = mean(Value))
    
    
    df <- bind_cols(rna,prot)
    
    sp1 <- ggscatter(
        df,
        x = "mean_TPM...2",
        y = 'mean_TPM...4',
        add = "reg.line",
        add.params = list(color = "blue", fill = "lightgray"),
        conf.int = TRUE
    )  +
        scale_x_log10() +
        scale_y_log10()+
        labs(x = 'RNA',y = 'Protein') +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),method = method,size = 10) +
        theme_bw() +
        theme(axis.text.x = element_text(size =15),
              axis.text.y = element_text(size =10),
              axis.title.y = element_text(size =15)) 
    sp1
}


NB4_gene_exp <- readRDS('Rds/nb4_gene_exp.Rds') 
NB4_gene_exp_DT <- 
    datatable(
        NB4_gene_exp,
        rownames = FALSE,
        options = list(
        autoWidth = TRUE,
        lengthChange = TRUE,
        pageLength = 30,
        searchHighlight = TRUE
    )) %>%
    formatRound(columns = str_subset(colnames(NB4_gene_exp), 'A\\d'), digits = 3)


NB4_prot_exp <- readRDS('Rds/nb4_protein_exp.Rds') 
NB4_prot_exp_DT <- 
    datatable(
        NB4_prot_exp,
        rownames = FALSE,
        extensions = 'FixedColumns',
        options = list(
            autoWidth = TRUE,
            pageLength = 30,
            scrollX = TRUE,
            searchHighlight = TRUE,
            fixedColumns = list(leftColumns =1)
        )) %>%
    formatRound(columns = str_subset(colnames(NB4_prot_exp), 'A\\d'), digits = 3)




header <- dashboardHeader(
    title = "NB4 viewer"
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("RNA expression", tabName = "RNA"),
        menuItem("Protein expressiobn", tabName = "Protein"),
        menuItem("Expression profie", tabName = "exp"),
        div(style="text-align:center")
    )
    
)

body <- dashboardBody(
    uiOutput("myTheme"),
    ##change the CSS format of sidebar 
    tags$head( 
        tags$style(HTML(".main-sidebar { font-size: 20px };")) #change the font size to 20
    ),
    
    tabItems(
        # RNA tab content
        tabItem(tabName = "RNA",
                div(
                    DT::dataTableOutput("rna_exp"),
                    style = "font-size:120%"
                )
        ),
        # Protein tab content
        tabItem(tabName = "Protein",
                div(
                    DT::dataTableOutput("prot_exp"),
                    style = "font-size:100%"
                )
        ),
        tabItem(tabName = "exp",
                box(
                    title = "Inputs gene name",
                    status = "primary", solidHeader = TRUE,
                    width = 12,
             
                    textInput('input_gene',NULL,value = 'MYC')
                ),
                fluidRow(
                    box(
                        title = "RNA expression",
                        status = "warning",
                        height = "800px", 
                        div(uiOutput("RNA_profile"),style = "margin-left: auto; margin-right: auto;")
                    ) ,
                    box(
                        title = "Protein expression",
                        status = "success",
                        height = "800px", 
                        div(uiOutput("protein_profile"),style = "margin-left: auto; margin-right: auto;")
                    ) 
                ),
                fluidRow(
                    box(
                        title = "Pearson correaltion between RNA and Protein",
                        status = "warning",
                        height = "800px", 
                        div(uiOutput("pearson"),style = "margin-left: auto; margin-right: auto;")
                    ) ,
                    box(
                        title = "Spearman correaltion between RNA and Protein",
                        status = "success",
                        height = "800px", 
                        div(uiOutput("sparman"),style = "margin-left: auto; margin-right: auto;")
                    ) 
                ),
            
        )
    )    
)

ui <-  dashboardPage(header, sidebar, body)

server <- function(input, output) {
    output$rna_exp <- renderDT(
        NB4_gene_exp_DT,
        server = TRUE,
        options = list(
            lengthChange = TRUE,
            pageLength = 30,
            searchHighlight = TRUE)
    )
    output$prot_exp <- renderDT(
        NB4_prot_exp_DT,
        server = TRUE,
        extensions = 'FixedColumns',
        options = list(
            dom = 'Bftip',
            pageLength = 30,
            searchHighlight = TRUE)
    )
    
    output$RNA_profile <- renderUI({
        if (toupper(input$input_gene) %in% NB4_gene_exp$Gene) {
            plotOutput("show_RNA_plot",height = "600px")
        } else {
            imageOutput("show_not_found")
        }
    })
    
    output$show_not_found <- renderImage({
        list(
            src = file.path("gene_not_found.png"),
            contentType = "image/png"
        )
    }, deleteFile = FALSE)
    
    output$protein_profile <- renderUI({
        if (toupper(input$input_gene) %in% NB4_prot_exp$Protein) {
            plotOutput("show_prot_plot",height = "600px")
        } else {
            imageOutput("show_not_found")
        }
    })
    
    output$show_RNA_plot <- renderPlot({
        plot_rna_exp(toupper(input$input_gene))
    }
    )
    output$show_prot_plot <- renderPlot({
        plot_prot_exp(toupper(input$input_gene))
    }
    )
    
    
    
    output$pearson <- renderUI({
        if (toupper(input$input_gene) %in% intersect(NB4_gene_exp$Gene,NB4_prot_exp$Protein)) {
            plotOutput("show_pearson_cor",height = "600px")
        } else {
            imageOutput("unpaired")
        }
    })
    output$unpaired <- renderImage({
        list(
            src = file.path("insuffect.png"),
            contentType = "image/png"
        )
    }, deleteFile = FALSE)

    
    output$show_pearson_cor <- renderPlot({
        plot_cor(toupper(input$input_gene),method = 'pearson')
    }
    )
    
    
    
    output$sparman <- renderUI({
        if (toupper(input$input_gene) %in% intersect(NB4_gene_exp$Gene,NB4_prot_exp$Protein)) {
            plotOutput("show_sparman_cor",height = "600px")
        } else {
            imageOutput("unpaired")
        }
    })
    
    output$show_sparman_cor <- renderPlot({
        plot_cor(toupper(input$input_gene),method = 'spearman')
    }
    )
    
    
    
}

shinyApp(ui, server)
