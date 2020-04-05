
library(shiny)
library(dplyr)
library(data.table)
library(plotly)
library(mvmeta)
library(shinythemes)

#options(error= recover)

load("dashData.Rdata")

##### Define server logic ##### 
server <- function(input, output) {
  
  ###### THE STATE AND YEAR SPECIFIC STUFF #####
  
  ##### set up menus ##### 
  output$yearmenu <- renderUI({
    selectInput(
      "year",
      "Select year",
      years,
      years[length(years)]
      
    )
  })
  
  output$popmenu <- renderUI({
    selectInput(
      "pop",
      "Select population",
      c("Total Population" = "_tot",
        "Males" = "_m",
        "Females" = "_f")
      
    )
  })
  
  populationname <- reactive({ #creates a reactive object that changes with input
    if (input$pop == "_tot") {
      popname <- "Total population"
    }
    if (input$pop == "_m") {
      popname <- "Male population"
    }
    if (input$pop == "_f") {
      popname <- "Female population"
    }
    popname
  })
  
  output$crimemenu <- renderUI({
    selectInput(
      "crimechoice",
      "Select crime",
      menus$Label
    )
  })
  
  output$statemenu <- renderUI({
    selectInput(
      "state",
      "Select state",
      states,
      "All States"
    )
  })
  
  output$customage <- renderUI({
    numericInput (
      "cusage",
      "Enter age",
      value = 15,
      min = 10,
      max = 90
    )
  })
  
  output$showtabbutton <- renderUI({
    radioButtons (
      "showtab",
      "Show tables below graphs?",
      c("No","Yes")
    )
  })
  
  output$uselnagebutton <- renderUI({
    radioButtons (
      "uselnage",
      "Use ln(age) in regressions?",
      c("No","Yes"),
      "Yes"
    )
  })
  
  output$oldagebutton <- renderUI({
    radioButtons (
      "oldage",
      "Use age 49-99 in regressions?",
      c("No","Yes"),
      "Yes"
    )
  })
  
  output$polymenu <- renderUI({
    if (nrow(crimedata())>4) {
      numericInput(
        "poly",
        "Number of terms in polynomial regression (k)",
        min = 1,
        max = 5,
        value = 3
      )
    }
  })
  
  ##### create data for analysis ##### 
  
  crime <- reactive({
    crime <- c(menus[which(menus$Label == input$crimechoice),"Var"])
  })
  
  data_vars <- reactive({
    data_vars <- c("agemid", 
                   "lnagemid",
                   "agecat",
                   paste0("r_",crime(),input$pop),
                   paste0("l_",crime(),input$pop)
    )
    return(data_vars)
  })
  
  crimedata <- reactive({
    req(input$crimechoice, input$pop, input$state, input$year, input$oldage)
    newdata <- dashdata %>% arrange(agemid) %>%
      filter(state==input$state) %>%
      filter(year==input$year) %>%
      select(data_vars()) %>%
      filter(is.infinite(!!sym(paste0("l_",crime(),input$pop))) == FALSE)
    
    if (input$oldage == "No") {
      newdata <- newdata %>% filter(agecat != "49-99")
    }
    as.data.frame(newdata)
  })
  
  ##### write out model #####
  
  mathtext <- reactive({
    form <- "$$\\Lambda\\left(\\frac{\\mbox{Arrests}_{age}}{\\mbox{Population}_{age}}\\right)=\\beta_0"
    if (input$uselnage == "Yes") {
      form <- paste0(form,"+\\beta_1 \\mbox{ln}(age)")
      if (input$poly > 1) {
        for (i in 2:input$poly) {
          form <- paste0(form,"+\\beta_",i,"\\mbox{ln}^",i," (age)")
        }
      }
    }
    else {
      form <- paste0(form,"+\\beta_1age")
      if (input$poly > 1) {
        for (i in 2:input$poly) {
          form <- paste0(form,"+\\beta_",i," age^",i)
        }
      }
    }
    form <- paste0(form,"+\\epsilon_{age}$$")
    return(form)
  })
  
  output$mathtextoutput <- renderUI({
    req(input$crimechoice, input$pop, input$state, input$year, input$oldage,input$showtab,input$uselnage)
    withMathJax(h5(mathtext()))
  })
  
  output$r2title <- renderText({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    paste("model r-square =",
          model_r_square())
  })
  
  ##### write out effect size #####
  
  effectsizemathext <- reactive({
    if (input$uselnage == "Yes") {
      form <- "$$d_{age} = \\frac{\\sqrt{3}}{\\pi} \\frac{1}{age} \\sum_k k \\beta_k \\mbox{ln}^{k-1}(age)$$"
    }
    else {
      form <- "$$d_{age} = \\frac{\\sqrt{3}}{\\pi}\\sum_k k \\beta_k age^{k-1}$$"
    }
    return(form)
  })
  
  output$effectsizemathtextoutput <- renderUI({
    req(input$crimechoice, input$pop, input$state, input$year, input$oldage,input$showtab,input$uselnage)
    withMathJax(h5(effectsizemathext()))
  })
  
  ##### create model (OLS) ##### 
  
  model_formula <- reactive({
    if (input$uselnage == "Yes") {
      form <- paste0("l_",crime(),input$pop,"~lnagemid")
      if (as.numeric(input$poly) > 1) {
        for (i in 2:as.numeric(input$poly)) {
          form <- paste0(form,"+I(lnagemid^",i,")")
        }
      }
    }
    else {
      form <- paste0("l_",crime(),input$pop,"~agemid")
      if (as.numeric(input$poly) > 1) {
        for (i in 2:as.numeric(input$poly)) {
          form <- paste0(form,"+I(agemid^",i,")")
        }
      }
    }
    return(form)
  })
  
  results_model <- reactive({ #estimate model
    lm(as.formula(model_formula()), data = crimedata())
  })
  
  #### Model fitted values ####
  
  model_fitted_rates <- reactive({
    1e+5*exp(results_model()$fitted.values)/
      (1+exp(results_model()$fitted.values))
  })
  
  
  ##### create effect sizes ##### 
  
  model_effect_sizes <- reactive({
    vals <- rep(0,nrow(crimedata()))
    for (i in as.numeric(input$poly):1) {
      if (is.na(results_model()$coefficients[i+1]) == FALSE) {
        if (input$uselnage == "Yes") {
          vals <- vals + 
            i*results_model()$coefficients[i+1]*log(crimedata()$agemid)^(i-1)/crimedata()$agemid
        }
        else {
          vals <- vals + 
            i*results_model()$coefficients[i+1]*crimedata()$agemid^(i-1)
        }
      }
    }
    vals <- sqrt(3)/pi*(vals)
    return(round(vals, digits = 3))
  })
  
  custom_effect_size <- reactive({
    vals <- 0
    for (i in as.numeric(input$poly):1) {
      if (is.na(results_model()$coefficients[i+1]) == FALSE) {
        if (input$uselnage == "Yes") {
          vals <- vals + 
            i*results_model()$coefficients[i+1]*log(input$cusage)^(i-1)/input$cusage
        }
        else {
          vals <- vals + 
            i*results_model()$coefficients[i+1]*input$cusage^(i-1)
        }
      }
    }
    vals <- sqrt(3)/pi*(vals)
    names(vals) <- paste("Effect size for age",input$cusage)
    return(round(vals, digits = 3))
  })
  
  
  #### Model Summaries #### 
  
  ##### create model plot (OLS) ##### 
  
  output$rawplot <- renderPlot({
    if (input$uselnage == "Yes") {
      x_var <- crimedata()$lnagemid
      x_lab <- "ln of age category midpoint"
    }
    else {
      x_var <- crimedata()$agemid
      x_lab <- "age category midpoint"
    }
    hats <- results_model()$fitted.values
    minval <- min(min(hats),min(crimedata()[,5]))
    maxval <- max(max(hats),max(crimedata()[,5]))
    plot(x_var,hats,type = "l",
         ylim = c(minval,maxval),
         ylab = "logit of rate*1e-5",
         xlab = x_lab
    )
    par(new=TRUE)
    plot(x_var,crimedata()[,5],
         ylim = c(minval,maxval),
         ylab = "logit of rate*1e-5",
         xlab = x_lab
    )
  })
  
  
  model_r_square <- reactive({
    as.character(
      round(
        summary(results_model())$r.squared
        ,
        digits = 3
      )
    )
  })
  
  output$modelsummary <- renderPrint({
    summary(results_model())
  })
  
  
  
  ##### create tables ##### 
  
  output$customeffectsize <- renderPrint({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    custom_effect_size()
  })
  
  rategraphdata <- reactive({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    cbind(crimedata(),model_fitted_rates())
  })
  
  output$rategraphdatatab <- renderTable({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    if (input$showtab == "Yes") {
      rategraphdata()[,c(-1,-2,-5)]
    }
    
  })
  
  effectsizegraphdata <- reactive({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    cbind(
      crimedata(),model_effect_sizes()
    )
  })
  
  output$effectsizegraphdataoutput <- renderTable({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    if (input$showtab == "Yes") {
      effectsizegraphdata()[,c(-1,-2,-4,-5)]
    }
    
  })
  
  
  output$alldatatab <- renderTable({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    newdata <- cbind(
      crimedata()[,c(-1,-2,-5)],model_fitted_rates(),model_effect_sizes()
    )
    colnames(newdata) <- c("Age category",
                           paste("Rate of",input$crime,"per 100k"),
                           paste("Predicted rate of",input$crime,"per 100k"),
                           "Age-based d effect size")
    newdata
  })
  
  ##### create rate plot ##### 
  
  rateplot <- reactive({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    graphvars <- colnames(rategraphdata())
    p <- plot_ly(rategraphdata(),
                 x = ~get(graphvars[1]),
                 y = ~get(graphvars[4]),
                 name = 'Actual rate',
                 type = 'scatter',
                 mode = 'lines',
                 line = list(color = "black"),
                 text = ~paste("Ages", get(graphvars[3]),
                               "\nActual rate = ", round(get(graphvars[4]), digits = 1)),
                 hoverinfo = "text") %>%
      add_trace(y = ~get(graphvars[6]),
                name = 'Model-based prediction of rate',
                mode = 'lines',
                line = list(color = "blue"),
                text = ~paste("Ages", get(graphvars[3]),
                              "\nPredicted rate = ", round(get(graphvars[6]), digits = 1)),
                hoverinfo = "text") %>%
      layout(xaxis = list(title = "Age"),
             yaxis = list (title = "Rate per 100,000"),
             legend = list(orientation = 'h',y = -0.25),
             plot_bgcolor = "white",
             paper_bgcolor = "white") %>%
      config(displayModeBar = TRUE)
    return(p)
  })
  
  output$rateplot <- renderPlotly({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    rateplot()
  })
  
  ##### create es plot ##### 
  
  effectsizeplot <- reactive({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$showtab)
    req(input$uselnage)
    graphvars <- colnames(effectsizegraphdata())
    p <- plot_ly(effectsizegraphdata(),
                 x = ~get(graphvars[1]),
                 y = ~get(graphvars[6]),
                 text = ~paste("Ages", get(graphvars[3]),
                               "\nd = ", get(graphvars[6])),
                 hoverinfo = "text",
                 name = 'Model-based prediction',
                 mode = 'markers+lines',
                 line = list(color = "blue")) %>%
      layout(xaxis = list(title = "Age"),
             yaxis = list (title = "d-type effect size")) %>%
      config(displayModeBar = TRUE)
    return(p)
  })
  
  output$ratestitle <- renderText({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    paste("Age based rates of",
          input$crime,
          "in \n",
          input$year,
          "for",
          input$state,
          "(",populationname(),")")
  })
  
  output$effectsizeplot <- renderPlotly({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    effectsizeplot()
  })
  
  output$effectsizetitle <- renderText({
    req(input$crimechoice)
    req(input$pop)
    req(input$state)
    req(input$year)
    req(input$poly)
    req(input$showtab)
    req(input$uselnage)
    paste("Age based effects sizes for",
          input$crime,
          "in \n",
          input$year,
          "for",
          input$state,
          "(",populationname(),")")
  })
  
  ###### THE META ANALYSIS #####
  
  output$meta_popmenu <- renderUI({
    selectInput(
      "meta_pop",
      "Select population",
      c("Total Population" = "_tot",
        "Males" = "_m",
        "Females" = "_f")
      
    )
  })
  
  output$meta_polymenu <- renderUI({
    numericInput(
      "meta_poly",
      "Number of terms in polynomial regression (k)",
      min = 1,
      max = 3,
      value = 3
    )
  })
  
  meta_populationname <- reactive({
    if (input$meta_pop == "_tot") {
      meta_popname <- "Total population"
    }
    if (input$meta_pop == "_m") {
      meta_popname <- "Male population"
    }
    if (input$meta_pop == "_f") {
      meta_popname <- "Female population"
    }
    meta_popname
  })
  
  output$meta_crimemenu <- renderUI({
    selectInput(
      "meta_crimechoice",
      "Select crime",
      menus$Label
    )
  })
  
  output$meta_oldagebutton <- renderUI({
    radioButtons (
      "meta_oldage",
      "Use age 49-99 in meta analysis?",
      c("No","Yes"),
      "Yes"
    )
  })
  
  #### Make Data ####
  
  meta_crime <- reactive({
    meta_crime <- c(menus[which(menus$Label == input$meta_crimechoice),"Var"])
  })
  
  meta_data_vars <- reactive({
    meta_data_vars <- c("agecat",
                        "agemid",
                        "lnagemid",
                        "state",
                        "year",
                        "pop_tot",
                        paste0(meta_crime(),input$meta_pop),
                        paste0("r_",meta_crime(),input$meta_pop),
                        paste0("l_",meta_crime(),input$meta_pop))
    return(meta_data_vars)
  })
  
  meta_crimedata <- reactive({
    req(input$meta_crimechoice, input$meta_pop, input$meta_oldage)
    newdata <- dashdata %>%
      filter(state != "All States") %>%
      filter(year != "Pooled 2000-2015") %>%
      filter(is.infinite(!!sym(paste0("l_",meta_crime(),input$meta_pop))) == FALSE) %>%
      select(meta_data_vars()) %>%
      mutate(age = agemid) %>%
      mutate(lnage = log(age))
    if (input$meta_oldage == "No") {
      newdata <- newdata %>% filter(agecat != "49-99")
    }
    as.data.frame(newdata)
  })
  
  list_of_states <- reactive({
    rownames(table(meta_crimedata()$state))
  })
  
  list_of_years <- reactive({
    rownames(table(meta_crimedata()$year))
  })
  
  agevals <- reactive({
    as.numeric(rownames(table(meta_crimedata()$age)))
  })
  
  meta_form <- reactive({
    form <- "cbind("
    for (p in 0:input$meta_poly) {
      form <- paste0(form,"b",p)
      if (p < input$meta_poly) {
        form <- paste0(form,",")
      }
    }
    as.formula(paste0(form,")~-1+factor(year)")) #note that there is no intercept
  })
  
  meta_analysis_model <- eventReactive(input$run_meta, {
    n <- length(list_of_years())*length(list_of_states())+1
    withProgress(message = "Running Meta Analysis", {
      b <- c()
      V <- list()
      for (y in list_of_years()) {
        for (s in 1:length(list_of_states())) {
          model_data <- meta_crimedata() %>%
            filter(state == list_of_states()[s]) %>%
            filter(year == y)
          if (nrow(model_data) == length(agevals())) {
            model <- lm(as.formula(paste0(paste0("l_",meta_crime(),input$meta_pop),
                                          "~poly(lnage, poly = input$meta_poly, raw = TRUE)")), 
                        data = model_data)
            coefs <- data.frame(t(coef(model)))
            colnames(coefs) <- paste0("b",0:input$meta_poly)
            coefs$year <- y
            coefs$state <- list_of_states()[s]
            b <- rbind(b, coefs)
            V[[paste0(s,y)]] <- vcov(model)
            incProgress(1/n, detail = paste("running OLS on ", list_of_states()[s], "-", y))
          }
        }
      }
      
      incProgress(1/n, detail = paste("running meta regression"))
      assign("meta_poly", input$meta_poly, envir = .GlobalEnv) #need to add poly to global env for some fucked up reason
      meta_mod <- mvmeta(as.formula(paste("cbind(",paste(paste0("b",0:meta_poly), collapse = ","),")~-1+factor(year)")),V, 
                         data = b)
      
    })
    return(meta_mod)
  })
  
  meta_coefs <- reactive({
    req(input$run_meta)
    coef(meta_analysis_model())
  })
  
  meta_coef_names <- reactive({
    req(input$run_meta)
    names(meta_coefs())
  })
  
  meta_V <- reactive({
    req(input$run_meta)
    vcov(meta_analysis_model())
  })
  
  meta_results <- reactive({
    poly <- 3
    results <- c()
    for (y in list_of_years()) {
      meta_logit <- rep(0, length(agevals()))
      for (p in 0:poly) {
        meta_logit <- meta_logit + 
          meta_coefs()[paste0("b",p,".factor(year)",y)]*log(agevals())^p
      }
      
      meta_effect_sizes <- c()
      meta_effect_sizes_se <- c()
      
      for (a in agevals()) {
        w <- rep(0,length(meta_coefs()))
        for (p in 0:poly) {
          w[which(meta_coef_names() == paste0("b",p,".factor(year)",y))] <- 
            sqrt(3)/pi*p*log(a)^(p-1)/a
        }
        meta_effect_sizes <- c(meta_effect_sizes, meta_coefs()%*%w)
        meta_effect_sizes_se <- c(meta_effect_sizes_se, sqrt(w%*%meta_V()%*%w))
      }
      
      results <- rbind(results,
                       cbind(agemid = agevals(),
                             meta_logit, 
                             meta_rate = 1e+5*exp(meta_logit)/(1+exp(meta_logit)),
                             meta_effect_sizes,
                             meta_effect_sizes_se,
                             year = as.numeric(y))
      )
    }
    
    check_data <- meta_crimedata() %>%
      group_by(age, agemid, year) %>%
      summarise(rate = mean(!!sym(paste0("r_",meta_crime(),input$meta_pop))),
                logit = mean(!!sym(paste0("l_",meta_crime(),input$meta_pop)))) 
    
    results <- merge(check_data, results)
    
    return(results)
  })
  
  output$meta_model_output <- renderPrint({
    req(input$run_meta)
    summary(meta_analysis_model())
  })
  
  meta_rate_plot <- reactive({
    r <- plot_ly(meta_results(),
                 y = ~agemid, 
                 x = ~as.numeric(year), 
                 z = ~rate, 
                 color = ~year,
                 type = 'scatter3d', mode = 'lines') 
    
    mr <- plot_ly(meta_results(),
                  y = ~agemid, 
                  x = ~as.numeric(year), 
                  z = ~meta_rate, 
                  line = list(dash = "dash"),
                  color = ~year,
                  type = 'scatter3d', mode = 'lines') 
    
    s1 <- subplot(r,mr, nrows = 1) %>% 
      layout(showlegend = FALSE,
             scene = list(
               xaxis = list(title = "Year"),
               yaxis = list(title = "Age"),
               zaxis = list(title = "Rate * 1e+5")
             ))
    return(s1)
  })
  
  output$render_meta_rate_plot <- renderPlotly({
    req(input$run_meta)
    meta_rate_plot()
  })
  
  meta_logit_plot <- reactive({
    l <- plot_ly(meta_results(),
                 y = ~agemid, 
                 x = ~as.numeric(year), 
                 z = ~logit, 
                 color = ~year,
                 type = 'scatter3d', mode = 'lines') 
    
    ml <- plot_ly(meta_results(),
                  y = ~agemid, 
                  x = ~as.numeric(year), 
                  z = ~meta_logit, 
                  line = list(dash = "dash"),
                  color = ~year,
                  type = 'scatter3d', mode = 'lines') 
    
    s1 <- subplot(l,ml, nrows = 1) %>% 
      layout(showlegend = FALSE,
             scene = list(
               xaxis = list(title = "Year"),
               yaxis = list(title = "Age"),
               zaxis = list(title = "Logit of rate")
             )) 
    return(s1)
  })
  
  output$render_meta_logit_plot <- renderPlotly({
    req(input$run_meta)
    meta_logit_plot()
  })
  
  meta_es_plot <- reactive({
    req(input$run_meta)
    plot_ly(meta_results(),
            y = ~agemid, 
            x = ~as.numeric(year), 
            z = ~meta_effect_sizes, 
            color = ~year,
            type = 'scatter3d', mode = 'lines') %>% 
      layout(showlegend = FALSE,
             scene = list(
               xaxis = list(title = "Year"),
               yaxis = list(title = "Age"),
               zaxis = list(title = "Effect size")
             ))
  })
  
  output$render_meta_es_plot <- renderPlotly({
    req(input$run_meta)
    meta_es_plot()
  })
  
  
  meta_effect_size_table <- reactive({
    newtable <- meta_results()
    newtable$meta_effect_sizes_test <- 
      newtable$meta_effect_sizes/newtable$meta_effect_sizes_se
    
    newtable$p <- 2*(1-pnorm(abs(newtable$meta_effect_sizes_test)))
    
    newtable$star <- ""
    
    newtable$star[which(newtable$p < .05)] <- "*"
    newtable$star[which(newtable$p < .01)] <- "**"
    newtable$star[which(newtable$p < .001)] <- "***"
    
    newtable$value <- paste(as.character(
      round(
        newtable$meta_effect_sizes, digits = 3
      )
    ),
    newtable$star
    )
    
    newtable <- data.table(newtable) %>% 
      dcast(age ~ year)
    
    return(newtable)
  })
  
  output$render_meta_effect_size_table <- renderTable({
    meta_effect_size_table()
  })
  
}

##### Define UI for application ##### 

ui <- navbarPage("Effect sizes for year-over-year changes in arrest rates \n BETA",
                 tabPanel("Analysis of specific years and states",
                          sidebarLayout(
                            sidebarPanel(
                              h3("Select data"),
                              uiOutput("popmenu"),
                              uiOutput("yearmenu"),
                              uiOutput("crimemenu"),
                              uiOutput("statemenu"),
                              h3("Model specification"),
                              uiOutput("polymenu"),
                              uiOutput("uselnagebutton"),
                              uiOutput("oldagebutton")
                            ),
                            mainPanel(
                              tabsetPanel(
                                type = "tabs",
                                tabPanel("Results",
                                         fluidRow(
                                           h2("Model-based effect size calculator"),
                                           uiOutput("customage"),
                                           verbatimTextOutput("customeffectsize")
                                         ),
                                         fluidRow(
                                           column(6,
                                                  h2("Model visualizations"),
                                                  h5(textOutput("r2title"))
                                           ),
                                           column(6,
                                                  uiOutput("showtabbutton")
                                           )
                                         ),
                                         fluidRow(
                                           column(6,
                                                  h4(textOutput("ratestitle")),
                                                  plotlyOutput("rateplot"), 
                                                  tableOutput("rategraphdatatab")
                                           ),    
                                           
                                           column(6,
                                                  h4(textOutput("effectsizetitle")),
                                                  plotlyOutput("effectsizeplot"),
                                                  tableOutput("effectsizegraphdataoutput")
                                           )
                                         )
                                ),
                                tabPanel("Model summary", 
                                         fluidRow(
                                           h4("Model estimated:"),
                                           uiOutput("mathtextoutput"),
                                           h4("Effect size estimated with:"),
                                           uiOutput("effectsizemathtextoutput"),
                                           h6("where k is the number of terms in the polynomial")
                                         ),
                                         fluidRow(
                                           column(6, 
                                                  h3("Statistical output:"),
                                                  verbatimTextOutput("modelsummary")
                                           ),
                                           column(6,
                                                  h3("Data plot:"),
                                                  plotOutput("rawplot"),
                                                  h5("points are empirical values and line is the fitted values")
                                           )
                                         )
                                ),
                                tabPanel("Data", 
                                         tableOutput("alldatatab")
                                )
                                
                              )
                            )
                          )
                 ),
                 tabPanel(
                   "Meta analysis across years and states",
                   sidebarLayout(
                     sidebarPanel(
                       h3("Select data"),
                       uiOutput("meta_popmenu"),
                       uiOutput("meta_crimemenu"),
                       uiOutput("meta_oldagebutton"),
                       uiOutput("meta_polymenu")
                     ),
                     mainPanel(
                       actionButton("run_meta", "Run meta analysis with selections"),
                       tabsetPanel(
                         type = "tabs",
                         tabPanel("Table of effect sizes",
                                  tableOutput("render_meta_effect_size_table")
                         ),
                         tabPanel("Effect sizes graph",
                                  plotlyOutput("render_meta_es_plot",height = 600, width = 600)
                         ),
                         tabPanel("Rate validation",
                                  plotlyOutput("render_meta_rate_plot",height = 600, width = 600),
                                  h4("solid line is average across states"),
                                  h4("dashed line is meta regression prediction")
                         ),
                         tabPanel("Logit validation",
                                  plotlyOutput("render_meta_logit_plot",height = 600, width = 600),
                                  h4("solid line is average across states"),
                                  h4("dashed line is meta regression prediction")
                         ),
                         tabPanel("Statistial output",
                                  verbatimTextOutput("meta_model_output")
                         )
                       )
                       
                     )
                   )
                 ),
                 tabPanel("Data sources",
                          h5("Kaplan, Jacob. Uniform Crime Reporting (UCR) Program Data: Arrests by Age, Sex, and Race, 1974-2016: ucr_arrests_yearly_1974_2016_dta.zip. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2018-12-29. https://doi.org/10.3886/E102263V7-10643"),
                          h5("Persistent URL:  http://doi.org/10.3886/E102263V7-10643"),
                          h4("and"),
                          h5("Steven Ruggles, Sarah Flood, Ronald Goeken, Josiah Grover, Erin Meyer, Jose Pacas and Matthew Sobek. IPUMS USA: Version 9.0 [dataset]. Minneapolis, MN: IPUMS, 2019. https://doi.org/10.18128/D010.V9.0")
                 )
                 
)

# Run the application 
shinyApp(ui = ui, server = server)
