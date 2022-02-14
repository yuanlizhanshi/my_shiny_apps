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
load('KWMT2NCBI.rdata')
load('CN_3G2NCBI.rdata')
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
ID_conversion <- function(df,TOtype){
  id_vec <- df[[1]]
  switch(
    TOtype,
    'KWMT2NCBI' = map_chr(id_vec,KWMT2NCBI_trans),
    'NCBI2KWMT' = map_chr(id_vec,NCBI2KWMT_trans),
    'CN_3G2NCBI' = map_chr(id_vec,CN_3G2NCBI_trans),
    'NCBI2CN_3G' = map_chr(id_vec,NCBI2CN_3G_trans)
  )
}

all_choices <- c('KWMT2NCBI','NCBI2KWMT','CN_3G2NCBI','NCBI2CN_3G')

getcolnames <- function(choices){
  case_when(choices == 'KWMT2NCBI' ~ c('Kaikobase_id','NCBI_id'),
            choices == 'NCBI2KWMT' ~ c('NCBI_id','Kaikobase_id'),
            choices == 'CN_3G2NCBI' ~ c('Silkdb3.0_id','NCBI_id'),
            choices == 'NCBI2CN_3G' ~ c('NCBI_id','Silkdb3.0_id'),
  )
}

ui <- fluidPage(
  fluidRow(
    fileInput("upload", NULL, buttonLabel = "Upload...", multiple = FALSE,accept = c(".xls", ".xlsx"))
  ),
  fluidRow(
    selectInput('convert_type','Select the conversion type',choices = all_choices)
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
)

server <- function(input, output, session) {
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
    convert_df <- data() %>% mutate(ID_conversion(data(),input$convert_type)) %>%
      data.table::setnames(getcolnames(input$convert_type))
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
}
shinyApp(ui = ui, server = server)

