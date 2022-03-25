library(shiny)
library(reticulate)
source_python('scihub.py')
download_paper <- function(papers){
  # run_splitstr.R
  command = "python"
  # Note the single + double quotes in the string (needed if paths have spaces)
  path2script='"scihub.py"'
  # Build up args in a vector
  pattern = "-d"
  args = c(pattern, papers)
  # Add path to script as first arg
  allArgs = c(path2script, args)
  output = system2(command, args=allArgs, stdout=TRUE,stderr = T)
  return(output)
}

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(0%);
             }"
      )
    )
  ),
  div(
    fluidRow(
      column(10,
             textInput('input_link','Downloads a paper from sci-hub given an indentifier (DOI, PMID, URL)',width = 1200)
      )
      
    ),style = "font-size:200%"
  ),
  actionButton("Search", "Search"),
  downloadButton("download_test", "Download paper")
)

server <- function(input, output, session) {
  observeEvent(input$Search,{
    download_res <- download_paper(input$input_link)
    if (any(grepl('Successfully',download_res))) {
      download_res <- 'Successfully parsed file '
      file.rename(list.files()[grep('*.pdf',list.files())],'output')
      file.remove(list.files()[grep('*.pdf',list.files())])
    }else{
      download_res <- 'Parsing file failed'
    }
    id <- showNotification(download_res,type = "message")
  })
  output$download <- downloadHandler(
    filename = function(){
      paste0(input$input_link, ".pdf")
    },
    content = function(file){
      file.copy("output", file)
    }
  )
}
shinyApp(ui = ui, server = server)
