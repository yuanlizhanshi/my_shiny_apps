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
library(AnnotationHub)
library(GO.db)

load('KWMT2NCBI.rdata')
load('CN_3G2NCBI.rdata')
load('all_ids.rdata')
load('KEGG_info.rdata')
load('hub.rdata')
load('domain_info.rdata')
###load database
hub <- AnnotationHub::AnnotationHub()
bmor_orgdb <- hub[['AH97101']]
GO <- as.list(GOTERM)
find_gene <- function(x){
  if(str_detect(x,'KWMT')){
    NCBI_ID_res <- KWMT2NCBI[str_which(KWMT2NCBI$KaikobaseID,x),]$NCBI_ID
    if (length(NCBI_ID_res) > 0) {
      id_type <- 'KaikobaseID'
      NCBI_ID <- NCBI_ID_res
      names(NCBI_ID) <- rep(id_type,length(NCBI_ID))
    }
  }else if(str_detect(x,'BMgn')){
    NCBI_ID_res <- KWMT2NCBI[str_which(KWMT2NCBI$old_gene_id,x),]$NCBI_ID
    id_type <- 'SilkbaseID'
    NCBI_ID <- NCBI_ID_res
    names(NCBI_ID) <- rep(id_type,length(NCBI_ID))
  }else if(str_detect(x,'BMSK')){
    NCBI_ID_res <- CN_3G2NCBI[str_which(CN_3G2NCBI$silkDB3.0_ID,x),]$SYMBOL
    if (length(NCBI_ID_res) > 0) {
      id_type <- 'Silkdb3.0ID'
      NCBI_ID <- NCBI_ID_res
      names(NCBI_ID) <- rep(id_type,length(NCBI_ID))
    }
  }else if(length(names(all_ids[str_which(all_ids,fixed(x,ignore_case=TRUE))])) > 0){
    id_type <- unique(names(all_ids[str_which(all_ids,fixed(x,ignore_case=TRUE))]))
    if (length(id_type)>1 | 'alias' %in% id_type  ) {
      if(any(which(nchar(all_ids[str_which(all_ids,fixed(x,ignore_case=TRUE))]) == nchar(x)))) {
        id_type <- names(which(nchar(all_ids[str_which(all_ids,fixed(x,ignore_case=TRUE))]) == nchar(x)))
        if ('alias' %in% id_type) {
          if (length(CN_3G2NCBI[str_which(CN_3G2NCBI$Gene_name,regex(fixed(paste0('^',x,'$')),ignore_case=TRUE)),]$SYMBOL) != 0) {
            NCBI_ID <- CN_3G2NCBI[str_which(CN_3G2NCBI$Gene_name,regex(fixed(paste0('^',x,'$')),ignore_case=TRUE)),]$SYMBOL
            names(NCBI_ID) <- rep('alias',length(NCBI_ID))
          }else{
            NCBI_ID <- 'NA'
            names(NCBI_ID) <- 'Gene not found'
          }
        }else{
          res <- AnnotationDbi::select(bmor_orgdb,keys = x,columns  = 'SYMBOL',keytype = id_type)
          NCBI_ID <- res[names(res) =='SYMBOL'] %>% as.character()
          names(NCBI_ID) <- rep(id_type,length(NCBI_ID))
        }
      }else{
        NCBI_ID <- 'NA'
        names(NCBI_ID) <- 'Gene not found'
      }
    }else if('SYMBOL' %in% id_type){
      NCBI_ID <- unique(all_ids[str_which(all_ids,fixed(x,ignore_case=TRUE))])
      names(NCBI_ID) <- rep('NCBI_id',length(NCBI_ID))
    }
    else{
      res <- AnnotationDbi::select(bmor_orgdb,keys = x,columns  = 'SYMBOL',keytype = id_type)
      NCBI_ID <- res[names(res) =='SYMBOL'] %>% as.character()
      names(NCBI_ID) <- rep(id_type,length(NCBI_ID))
    }
  }
  else{
    NCBI_ID <- 'NA'
    names(NCBI_ID) <- 'Gene not found'
  }
  NCBI_ID <- NCBI_ID[!duplicated(NCBI_ID)]
  return(NCBI_ID)
}
add_info <- function(chr){
  if(any(is.na(chr))){
    chr <- '暂无，欢迎补充'
    return(chr)
  }else{
    return(chr)
  }
}
###load function


# Define UI for application that draws a histogram
ui <- fluidPage(
  helpText('输入任意你想查询的基因id,自动检测并支持多种形式，例如BMSK0000001,KWMTBOMO00001,BMgn002073,LOC119629070,trx,NM_001042449.1,GO:0000003,
    也可以输入文献的pubmed号查询是否有相关基因，例如10066809'),
  fluidRow(
    column(2,textInput('input_Id','Input_Id:',value = 'Trx')),
  ),
  fluidRow(
    verbatimTextOutput('input_Id_type'),
  ),
  fluidRow(
    verbatimTextOutput('NCBI_id')
  ),
  fluidRow(
    verbatimTextOutput('Silkdb3.0id')
  ),
  fluidRow(
    verbatimTextOutput('KaikobaseID')
  ),
  fluidRow(
    verbatimTextOutput('SilkbaseID')
  ),
  fluidRow(
    verbatimTextOutput('ENTREZID')
  ),
  fluidRow(
    verbatimTextOutput('Alias')
  ),
  fluidRow(
    verbatimTextOutput('GENENAME')
  ),
  fluidRow(
    verbatimTextOutput('REFSEQ')
  ),
  fluidRow(
    verbatimTextOutput('Gene_location')
  ),
  fluidRow(
    verbatimTextOutput('GO_infomation'),
    div(
      tableOutput('GO'),
        style = "font-size:80%"
      )
  ),
  fluidRow(
    verbatimTextOutput('KEGG_infomation'),
    div(
      tableOutput('KEGG'),
      style = "font-size:80%"
    )
  ),
  fluidRow(
    verbatimTextOutput('Protein_domain_infomation'),
    div(
      tableOutput('domian'),
      style = "font-size:80%"
    ),
    div(
      uiOutput("Protein_domain_note"),
      style = "font-size:70%"
    ),
  ),
  fluidRow(
    verbatimTextOutput('PMID'),
  ),
  fluidRow(
    verbatimTextOutput('des'),
  ),
  fluidRow(
    verbatimTextOutput('Summary'),
    tags$style(type='text/css', '#txt_out {white-space: pre-wrap;}')
  ),
  fluidRow(
      div(
        selectInput("homolog_all", "  Choose species", colnames(KWMT2NCBI)[20:7],selected = 'Drosophila_melanogaster'),
        style = "font-size:80%"
      )
    ),
  fluidRow(
    verbatimTextOutput('homolog'),
  ),
  helpText('若发现数据有错误或者缺失，请联系孔某'),
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  )
)

# Define server logic required to draw a histogram
server <- function(input, output){
  NCBI_ID <- reactive(find_gene(req(input$input_Id)))
  output$input_Id_type <- renderText({
    paste0('输入的id类型为:  ',names(NCBI_ID()))
  })
  output$NCBI_id <- renderText({
    paste0('NCBI_id:  ',NCBI_ID())
  }) 
  output$Silkdb3.0id <- renderText({
    Silkdb3.0_gene <- CN_3G2NCBI[str_which(CN_3G2NCBI$SYMBOL,NCBI_ID()),]$silkDB3.0_ID
    paste0('Silkdb3.0_id:  ',paste0(Silkdb3.0_gene,collapse = ';')) %>% unique()
})
  output$KaikobaseID <- renderText({
    Kaikobase_gene <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),]$KaikobaseID
    paste0('Kaikobase_id:  ',paste0(Kaikobase_gene,collapse = ';')) %>% unique()
  })
  output$SilkbaseID <- renderText({
    SilkbaseID_gene <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),]$old_gene_id
    paste0('Silkbase_id:  ',paste0(SilkbaseID_gene,collapse = ';')) %>% unique()
  })
  output$ENTREZID <- renderText({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'ENTREZID',keytype = 'SYMBOL')
    ENTREZ_ID <- tran_res[names(tran_res) ==  'ENTREZID' ][[1]] %>% unique()
    paste0('ENTREZID:  ',paste0(ENTREZ_ID,collapse = ';'))
  })
  output$KEGG_infomation <- renderText({
    c('KEGG_infomation:')
  })
  output$KEGG <- renderTable({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'ENTREZID',keytype = 'SYMBOL')
    ENTREZ_ID <- tran_res[names(tran_res) ==  'ENTREZID' ][[1]] %>% unique()
    map_dfr(ENTREZ_ID,function(x)keggdata_id[str_which(keggdata_id$ENTREZID,x),-2])
  })
  output$Protein_domain_infomation <- renderText({
    paste0('Protein_domain_infomation')
  })
  output$Protein_domain_note <- renderUI({
    div(tags$p("Note:For more domain infomation: Please entry",style = "display: inline;"),
        a("Pfam", href="http://pfam.xfam.org/",style = "display: inline;"),
        tags$p("and",style = "display: inline;"),
        a("Interpro", href="https://www.ebi.ac.uk/interpro/",style = "display: inline;"),
        tags$p(' ')
    )
  })
  output$domian <- renderTable({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'ENTREZID',keytype = 'SYMBOL')
    ENTREZ_ID <- tran_res[names(tran_res) ==  'ENTREZID' ][[1]] %>% unique()
    domain_info %>% filter(Gene %in% ENTREZ_ID) 
  })
  output$Alias <- renderText({
    Alias_gene <- CN_3G2NCBI[str_which(CN_3G2NCBI$SYMBOL,NCBI_ID()),]$Gene_name
    paste0('Alias:  ',paste0(Alias_gene,collapse = ';')) %>% unique()
  })
  output$GENENAME <- renderText({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'GENENAME',keytype = 'SYMBOL')
    GENENAME_ID <- tran_res[names(tran_res) ==  'GENENAME' ][[1]] %>% unique()
    paste0('Full name:  ',paste0(GENENAME_ID,collapse = ';'))
  })
  output$REFSEQ <- renderText({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'REFSEQ',keytype = 'SYMBOL')
    REFSEQ_ID <- tran_res[names(tran_res) ==  'REFSEQ' ][[1]] %>% unique()
    paste0('REFSEQ:  ',paste0(REFSEQ_ID,collapse = ';'))
  })
  output$Gene_location <- renderText({
    Gene_location <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),]$location
    paste0('Gene location:  ',paste0(Gene_location,collapse = ';')) %>% unique()
  })
  output$GO_infomation <- renderText({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'GO',keytype = 'SYMBOL')
    GO_ID <- tran_res[names(tran_res) ==  'GO' ][[1]] %>% unique() %>% add_info()
    paste0('GO_ID:  ',paste0(GO_ID,collapse = ';'))
  })
  output$GO <- renderTable({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'GO',keytype = 'SYMBOL')
    GO_ID <- tran_res[names(tran_res) ==  'GO' ][[1]] %>% unique()
    GO_full_des <- map_dfr(GO_ID,function(x){
      GO_res <- GO[[x]]
      GO_res_df <- data.frame(
        GOID = GO_res@GOID,
        Term = GO_res@Term,
        Ontology = GO_res@Ontology,
        Definition = GO_res@Definition
      )
      return(GO_res_df)
    })
  })
  output$PMID <- renderText({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'PMID',keytype = 'SYMBOL')
    PMID_ID <- tran_res[names(tran_res) ==  'PMID' ][[1]] %>% unique()
    paste0('PMID(Related paper):  ',paste0(PMID_ID,collapse = ';'))
  })
  output$des <- renderText({
    Kaikobase_gene_des <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),]$Description %>% unique() %>% add_info()
    gene_des_res <- paste0('Gene description:  ',paste0(Kaikobase_gene_des,collapse = ';')) %>% unique() 
  })  
  output$Summary <- renderText({
    Kaikobase_gene_sum <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),]$Summary %>% unique() %>% add_info()
    gene_sum_res <- paste0('Gene summary:  ',paste0(Kaikobase_gene_sum,collapse = ';')) %>% unique()  %>% str_wrap(width = 220) 
  })
  output$homolog <- renderText({
    Kaikobase_gene <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,NCBI_ID()),][[input$homolog_all]]
    paste0('Homolog:  ',paste0(Kaikobase_gene,collapse = ';'))
  }
  )
}

# Run the application

shinyApp(ui = ui, server = server)




