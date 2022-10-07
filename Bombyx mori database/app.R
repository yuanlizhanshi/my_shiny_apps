library(shiny)
library(tidyverse)
library(shinydashboard)
library(shinyWidgets)
##load database------
load('KWMT2NCBI.rdata')
load('CN_3G2NCBI.rdata')
load('all_ids.rdata')
load('KEGG_info.rdata')
load('domain_info.rdata')
load('GO_info.rdata')
load('PMID_title.rdata')


######ID search function------
bmor_orgdb <-  AnnotationDbi::loadDb(file = "bmor_orgdb.sqlite")
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
    chr <- 'Gene not found'
    return(chr)
  }else{
    return(chr)
  }
}

######ID transformation function-------
load('KWMT2NCBI.rdata')
load('CN_3G2NCBI.rdata')
load('NCBI_info_res.rdata')
KWMT2NCBI_trans <- function(x){
  temp_res <- KWMT2NCBI[str_which(KWMT2NCBI$KaikobaseID,x),]$NCBI_ID
  if (length(temp_res) >=1) {
    return(paste0(temp_res,collapse = ','))
  }
  else{
    return(NA)
  }
}
NCBI2KWMT_trans <- function(x){
  temp_res <- KWMT2NCBI[str_which(KWMT2NCBI$NCBI_ID,x),]$KaikobaseID
  if (length(temp_res) >=1) {
    return(paste0(temp_res,collapse = ','))
  }
  else{
    return(NA)
  }
}
CN_3G2NCBI_trans <- function(x){
  temp_res <- CN_3G2NCBI[str_which(CN_3G2NCBI$silkDB3.0_ID,x),]$SYMBOL
  if (length(temp_res) >=1) {
    return(paste0(temp_res,collapse = ','))
  }
  else{
    return(NA)
  }
}
NCBI2CN_3G_trans <- function(x){
  temp_res <- CN_3G2NCBI[str_which(CN_3G2NCBI$SYMBOL,x),]$silkDB3.0_ID
  if (length(temp_res) >=1) {
    return(paste0(temp_res,collapse = ','))
  }
  else{
    return(NA)
  }
}
SYMBOL2GENENAME <- function(x){
  temp_res <- NCBI_info[str_which(NCBI_info$SYMBOL,fixed(x)),]$GENENAME
  if (length(temp_res) >=1) {
    return(paste0(temp_res,collapse = ','))
  }
  else{
    return(NA)
  }
}

ID_conversion <- function(df,TOtype,session = session){
  id_vec <- df[[1]]
  switch(
    TOtype,
    'Kaikobase ID to NCBI ID' = map_chr(1:length(id_vec),function(x){
      updateProgressBar(
        session = session,
        id = "pb",
        value = x, total = length(id_vec),
        title = paste("Process", 'finished')
      )
      x = id_vec[x]
      return(KWMT2NCBI_trans(x))
    }),
    'NCBI ID tO Kaikobase ID' = map_chr(1:length(id_vec),function(x){
      updateProgressBar(
        session = session,
        id = "pb",
        value = x, total = length(id_vec),
        title = paste("Process", 'finished')
      )
      x = id_vec[x]
      return(NCBI2KWMT_trans(x))
    }),
    'Silkdb3.0 ID to NCBI ID' = map_chr(1:length(id_vec),function(x){
      updateProgressBar(
        session = session,
        id = "pb",
        value = x, total = length(id_vec),
        title = paste("Process", 'finished')
      )
      x = id_vec[x]
      return(CN_3G2NCBI_trans(x))
    }),
    'NCBI ID tO Silkdb3.0 ID' = map_chr(1:length(id_vec),function(x){
      updateProgressBar(
        session = session,
        id = "pb",
        value = x, total = length(id_vec),
        title = paste("Process", 'finished')
      )
      x = id_vec[x]
      return(NCBI2CN_3G_trans(x))
    })
  )
}
all_choices <- c('Kaikobase ID to NCBI ID','NCBI ID tO Kaikobase ID','Silkdb3.0 ID to NCBI ID','NCBI ID tO Silkdb3.0 ID')
getcolnames <- function(choices){
  case_when(choices == 'Kaikobase ID to NCBI ID' ~ c('Kaikobase_id','NCBI_id'),
            choices == 'NCBI ID to Kaikobase ID' ~ c('NCBI_id','Kaikobase_id'),
            choices == 'Silkdb3.0 ID to NCBI ID' ~ c('Silkdb3.0_id','NCBI_id'),
            choices == 'NCBI ID to Silkdb3.0 ID' ~ c('NCBI_id','Silkdb3.0_id'),
  )
}

######GO&KEGG enrichment function-------
KEGG_enrich <- function(x) {
  mapped.gene <- clusterProfiler::bitr(
    geneID = x,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = bmor_orgdb
  )
  enrich1 <- clusterProfiler::enrichKEGG(gene = mapped.gene$ENTREZID,
                                         keyType = "kegg",
                                         organism  = 'bmor',
                                         pvalueCutoff = 1,
                                         use_internal_data = T)
  enrich2 <- clusterProfiler::setReadable(enrich1, OrgDb = bmor_orgdb, keyType="ENTREZID")
  return(enrich2@result)
}
GO_enrich <- function(x){
  enrich <- clusterProfiler::enrichGO(gene = x,
                                      OrgDb = bmor_orgdb,
                                      keyType = "SYMBOL",
                                      ont = "ALL",
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
  return(enrich@result)
}
Enrich_conversion <- function(df,TOtype){
  id_vec <- df[[1]]
  switch(
    TOtype,
    'KEGG enrichment' = KEGG_enrich(id_vec),
    'GO enrichment' = GO_enrich(id_vec),
  )
}

######Web front-------
header <- dashboardHeader(
  title = "Silkworm database"
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("ID search", tabName = "ID_search"),
    menuItem("ID barch transform", tabName = "ID_transform"),
    menuItem("GO&KEGG enrichment", tabName = "ID_enrichment"),
    div(style="text-align:center")
  )
  
)

body <- dashboardBody(
  ##change the CSS format of sidebar 
  tags$head( 
    tags$style(HTML(".main-sidebar { font-size: 20px };")) #change the font size to 20
  ),
  
  tabItems(
    #### ID_search tab content------
    
    tabItem(tabName = "ID_search",
            helpText('Input gene which you want to search, several froms are supported, sush as BMSK0000001,KWMTBOMO00001,BMgn002073,LOC119629070,trx,NM_001042449.1 ...'),
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
              verbatimTextOutput('paper_info'),
              div(
                tableOutput('PMID'),
                style = "font-size:80%"
              ),
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
            helpText('If you have any question,please contact kongyunhui1@gmail.com'),
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
            )
    ),
    
    #### ID_transform tab content------
    tabItem(tabName = "ID_transform",
            fluidRow(
              fileInput("upload", NULL, buttonLabel = "Upload...", multiple = FALSE,accept = c(".xls", ".xlsx"))
            ),
            fluidRow(
              selectInput('convert_type','Select the conversion type',choices = all_choices)
            ),
            fluidRow(
              column(
                3,
                progressBar(
                  id = "pb",
                  value = 0,
                  total = 100,
                  title = "",
                  display_pct = TRUE
                )
              )
            ),
            fluidRow(
              'Preview the results'
            ),
            div(
              fluidRow(
                column(6,
                       dataTableOutput("preview")
                )
              ),
              style = "font-size:80%"
            ),
            downloadButton("download", "Download Conversion result")
    ),
    #### ID_enrichment tab content------
    tabItem(tabName = "ID_enrichment",
            fluidRow(
              helpText('Please input NCBI id')
            ),
            fluidRow(
              fileInput("enrich_upload", NULL, buttonLabel = "Upload...", multiple = FALSE,accept = c(".xls", ".xlsx"))
            ),
            fluidRow(
              selectInput('convert_type2','Select the conversion type',choices = c('KEGG enrichment','GO enrichment'))
            ),
            fluidRow(
              'Preview the results'
            ),
            div(
              fluidRow(
                column(6,
                       dataTableOutput("Enrichment_preview")
                )
              ),
              style = "font-size:80%"
            ),
            downloadButton("Enrichment_download", "Download enrichment result")     
            
    )
  )
)

ui <-  dashboardPage(header, sidebar, body,
                     skin = "green")

######Web back-------
server <- function(input, output,session = session){
  #############sever of ID search--------
  NCBI_ID <- reactive(find_gene(req(input$input_Id)))
  observeEvent(input$input_Id,{
    if (input$input_Id != 'Trx') {
      write.table(paste(Sys.time(),'\t',input$input_Id),file = 'log.txt',append = T,col.names = F,row.names = F,quote = F)
    }
  })
  output$input_Id_type <- renderText({
    paste0('Input id type is:  ',names(NCBI_ID()))
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
  output$paper_info <- renderText({
    paste0('Relevant paper:')
  })
  output$PMID <- renderTable({
    tran_res <- AnnotationDbi::select(bmor_orgdb,keys = NCBI_ID() ,columns  = 'PMID',keytype = 'SYMBOL')
    PMID_ID <- tran_res[names(tran_res) ==  'PMID' ][[1]] %>% unique()
    PMID_res %>% filter(PMID == PMID_ID)
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
  #############sever of transformation--------
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           xls = readxl::read_xls(input$upload$datapath),
           xlsx = readxl::read_xlsx(input$upload$datapath),
           validate("Invalid file; Please upload a .xls or .xlsx file")
    )
  })
  convert_res <- reactive({
    req(input$convert_type)
    convert_df <- data() %>% mutate(ID_conversion(data(),input$convert_type,session = session)) %>%
      data.table::setnames(getcolnames(input$convert_type)) %>%
      mutate(Description = map_chr(.data[['NCBI_id']],SYMBOL2GENENAME))
  })
  
  output$preview <- renderDataTable({
    convert_res()
  },
  options = list(pageLength = 20)
  )
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$convert_type, ".xls")
    },
    content = function(file) {
      write.table(convert_res(), file,sep = '\t',col.names = T,row.names = F,quote = F)
    }
  )
  #############sever of enrichment--------
  enrich_data <- reactive({
    req(input$enrich_upload)
    ext <- tools::file_ext(input$enrich_upload$name)
    switch(ext,
           xls = readxl::read_xls(input$enrich_upload$datapath),
           xlsx = readxl::read_xlsx(input$enrich_upload$datapath),
           validate("Invalid file; Please upload a .xls or .xlsx file")
    )
  })
  enrichres <- reactive({
    req(input$convert_type2)
    enrichres_df <- enrich_data() %>% Enrich_conversion(input$convert_type2)
  })
  
  output$Enrichment_preview <- renderDataTable({
    enrichres()
  },
  options = list(pageLength = 20)
  )
  output$Enrichment_download <- downloadHandler(
    filename = function() {
      paste0(input$convert_type2,' results', ".xls")
    },
    content = function(file) {
      write.table(enrichres(), file,sep = '\t',col.names = T,row.names = F,quote = F)
    }
  )
}

shinyApp(ui, server)
