library(shiny)
library(shinydashboard)
library(dplyr)
library(glue)
library(shinyauthr)
library(RSQLite)
library(DBI)
library(lubridate)
library(readxl)
library(writexl)
library(DT)

load_user <- data.table::fread('user.info')
cookie_expiry <- 7

user_base <- tibble(
  user = load_user$user,
  password = load_user$pass_word,
  password_hash = sapply(load_user$pass_word, sodium::password_store),
  permissions = load_user$permissions,
  name = load_user$user
)

get_sessions_from_db <- function(conn = db, expiry = cookie_expiry) {
  dbReadTable(conn, "sessions") %>%
    mutate(login_time = ymd_hms(login_time)) %>%
    as_tibble() %>%
    filter(login_time > now() - days(expiry))
}

# This function must accept two parameters: user and sessionid. It will be called whenever the user
# successfully logs in with a password.

add_session_to_db <- function(user, sessionid, conn = db) {
  tibble(user = user, sessionid = sessionid, login_time = as.character(now())) %>%
    dbWriteTable(conn, "sessions", ., append = TRUE)
}


db <- dbConnect(SQLite(), ":memory:")
dbCreateTable(db, "sessions", c(user = "TEXT", sessionid = "TEXT", login_time = "TEXT"))


#Load data from excel
df = read_excel("data.xlsx")
# UI部分
ui <- dashboardPage(
  dashboardHeader(title = span("Antibody repertory",
                               style = " font-size: 14px"),
                  tags$li(
                    class = "dropdown",
                    style = "padding: 8px;",
                    shinyauthr::logoutUI("logout")
                  ),
                  tags$li(
                    class = "dropdown",
                    tags$a(
                      icon("github"),
                      href = "https://github.com/yuanlizhanshi",
                      title = "Show source code in kongmou's github"
                    )
                  )
                  
  ),

  
  dashboardSidebar(
    sidebarMenu(
      menuItem("试剂库存", tabName = "reagents", icon = icon("flask")),
      menuItem("取用登记", tabName = "logs", icon = icon("list-alt"))
    )
  ),
  
  dashboardBody(
    shinyauthr::loginUI(
      "login", 
      cookie_expiry = cookie_expiry, 
      additional_ui = tagList(
        tags$p("Welcome to wenlab anitbody management platform", class = "text-center",style = "padding: 20px"),
      )
    ),
    uiOutput("testUI")

  )
)

server <- function(input, output,session) {
  
  credentials <- shinyauthr::loginServer(
    id = "login",
    data = user_base,
    user_col = user,
    pwd_col = password_hash,
    sodium_hashed = TRUE,
    cookie_logins = TRUE,
    sessionid_col = sessionid,
    cookie_getter = get_sessions_from_db,
    cookie_setter = add_session_to_db,
    log_out = reactive(logout_init())
  )
  logout_init <- shinyauthr::logoutServer(
    id = "logout",
    active = reactive(credentials()$user_auth)
  )
  observe({
    if (credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })
  user_info <- reactive({
    credentials()$info
  })
  output$welcome <- renderText({
    req(credentials()$user_auth)
    glue("Welcome {user_info()$name}")
  })

  output$testUI <- renderUI({
    req(credentials()$user_auth)
    if (user_info()$permissions == "admin") {
      tabItem(tabName = "reagents",
              fluidRow(
                box(width = 12, title = "库存情况", status = "primary", solidHeader = TRUE,
                    dataTableOutput("table1"))
              ),
              fluidRow(
                box(width = 12, title = "添加试剂", status = "primary", solidHeader = TRUE,
                    selectInput("name", "试剂名称：", choices = unique(df$regent)),
                    numericInput("quantity1", "库存数量：", value = 2),
                    actionButton("add_button", "添加试剂"))
              )
      )
    }else if (user_info()$permissions == "standard") {
      tabItem(tabName = "logs",
              fluidRow(
                box(width = 12, title = "库存情况", status = "primary", solidHeader = TRUE,
                    dataTableOutput("table2"))
              ),
              fluidRow(
                box(width = 12, title = "取用登记", status = "primary", solidHeader = TRUE,
                    selectInput("reagent_name", "试剂名称：", unique(df$regent)),
                    numericInput("quantity2", "取用数量：", value = 1),
                    actionButton("log_button", "登记"))
              )
      )
    }
    
  })
  
  # reactive variables
  rv <- reactiveValues(df1 = df)
  
  
  observeEvent(input$add_button, {
    rv$df1[rv$df1$regent == input$name, "num"] <- rv$df1[rv$df1$regent == input$name, "num"] + input$quantity1
    write_xlsx(rv$df1, "data.xlsx")
  })

  output$table1 <- renderDataTable({
    rv$df1
  },
  options = list(pageLength = 10)
  )
  output$table2 <- renderDataTable({
    rv$df1
  },
  options = list(pageLength = 10)
  )
  
  observeEvent(input$log_button, {
    if (input$quantity2 > 0) {
        rv$df1[rv$df1$regent == input$reagent_name, "num"] <- rv$df1[rv$df1$regent == input$reagent_name, "num"] - input$quantity2
        logs_df <- data.frame(regent = input$reagent_name, quantity = input$quantity2, name = user_info()$name,date = Sys.time(), stringsAsFactors = FALSE)
        write.table(logs_df, file = "logs.txt", sep = '\t', append = TRUE,col.names = F,row.names = F,quote = F)
        write_xlsx(rv$df1, "data.xlsx")
      
    }
  })
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui, server)