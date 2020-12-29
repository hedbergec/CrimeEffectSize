rm(list = ls(all.names = TRUE))

library(data.table)
library(plotly)
library(mvmeta)
library(plot3D)


load("dashData.Rdata")

# yvar <- "totaldrug"
# yvar <- "stolenprop"
# yvar <- "mtrvehtheft" 
# yvar <- "weapons"
# yvar <- "disordercond"
yvar <- "aggassault"

poly <- 3

data_vars <- c("agecat",
               "agemid",
               "lnagemid",
               "state",
               "year",
               "pop_tot",
               paste0(yvar,"_tot"),
               paste0("r_",yvar,"_tot"),
               paste0("l_",yvar,"_tot"))

input_data <- dashdata %>%
  filter(state != "All States") %>%
  filter(year != "Pooled 2000-2015") %>%
  #filter(agemid < 65) %>%
  #filter(year > 1989) %>%
  filter(is.infinite(!!sym(paste0("l_",yvar,"_tot"))) == FALSE) %>%
  select(data_vars) %>%
  mutate(age = agemid) %>%
  mutate(lnage = log(age))

list_of_states <- rownames(table(input_data$state))
list_of_years <- rownames(table(input_data$year))

agevals <- as.numeric(rownames(table(input_data$age)))

b <- c()
V <- list()
for (y in list_of_years) {
  for (s in 1:length(list_of_states)) {
    model_data <- input_data %>%
      filter(state == list_of_states[s]) %>%
      filter(year == y)
    if (nrow(model_data) == length(agevals)) {
      model <- lm(as.formula(paste0(paste0("l_",yvar,"_tot"),
                                    "~poly(lnage,poly, raw = TRUE)")), 
                  data = model_data)
      coefs <- data.frame(t(coef(model)))
      colnames(coefs) <- paste0("b",0:poly)
      coefs$year <- y
      coefs$state <- list_of_states[s]
      b <- rbind(b, coefs)
      V[[paste0(s,y)]] <- vcov(model)
    }
  }
}

form <- "cbind("

for (p in 0:poly) {
  form <- paste0(form,"b",p)
  if (p < poly) {
    form <- paste0(form,",")
  }
}

form <- paste0(form,")~-1+factor(year)") #note that there is no intercept
meta_mod <- mvmeta(as.formula(form),V, data = b)
meta_coefs <- coef(meta_mod)
meta_coef_names <- names(meta_coefs)
meta_V <- vcov(meta_mod)

results <- c()


for (y in list_of_years) {
  meta_logit <- rep(0, length(agevals))
  for (p in 0:poly) {
    meta_logit <- meta_logit + meta_coefs[paste0("b",p,".factor(year)",y)]*log(agevals)^p
  }
  
  meta_effect_sizes <- c()
  meta_effect_sizes_se <- c()
  
  for (a in agevals) {
    w <- rep(0,length(meta_coefs))
    for (p in 0:poly) {
      w[which(meta_coef_names == paste0("b",p,".factor(year)",y))] <- 
        sqrt(3)/pi*p*log(a)^(p-1)/a
    }
    meta_effect_sizes <- c(meta_effect_sizes, meta_coefs%*%w)
    meta_effect_sizes_se <- c(meta_effect_sizes_se, sqrt(w%*%meta_V%*%w))
  }
  
  results <- rbind(results,
                   cbind(agemid = agevals,
                         meta_logit, 
                         meta_rate = 1e+5*exp(meta_logit)/(1+exp(meta_logit)),
                         meta_effect_sizes,
                         meta_effect_sizes_se,
                         year = as.numeric(y))
                   )
}


check_data <- input_data %>%
  group_by(age, agemid, year) %>%
  summarise(rate = mean(!!sym(paste0("r_",yvar,"_tot"))),
            logit = mean(!!sym(paste0("l_",yvar,"_tot")))) 

check_data <- merge(check_data, results)

p1 <- plot_ly(check_data, 
              y =~rate, 
              x =~agemid, 
              color =~year,
              legendgroup = ~year,
              type = "scatter" , 
              mode = 'lines') %>% layout(showlegend = FALSE, title = yvar,
                                             yaxis = list(range = c(min(check_data$rate), max(check_data$rate))))
p2 <- plot_ly(check_data, 
              y =~meta_rate, 
              x =~agemid, 
              color =~year,
              legendgroup = ~year,
              type = "scatter" , 
              mode = 'lines')%>% layout(showlegend = FALSE, title = yvar,
                                            yaxis = list(range = c(min(check_data$rate), max(check_data$rate))))

s1 <- subplot(p1,p2, nrows = 1) 

p1 <- plot_ly(check_data, 
              y =~logit, 
              x =~agemid, 
              color =~year,
              legendgroup = ~year,
              type = "scatter" , 
              mode = 'lines') %>% layout(showlegend = FALSE, title = yvar,
                                         yaxis = list(range = c(min(check_data$logit), max(check_data$logit))))
p2 <- plot_ly(check_data, 
              y =~meta_logit, 
              x =~agemid, 
              color =~year,
              legendgroup = ~year,
              type = "scatter" , 
              mode = 'lines')%>% layout(showlegend = FALSE, title = yvar,
                                        yaxis = list(range = c(min(check_data$logit), max(check_data$logit))))

s2 <- subplot(p1,p2, nrows = 1) 

c1 <- subplot(s1,s2, nrows = 2) 

e1 <- plot_ly(check_data, 
              y =~meta_effect_sizes, 
              x =~agemid, 
              color =~year,
              legendgroup = ~year,
              type = "scatter" , 
              mode = 'lines')%>% 
  add_ribbons(ymin = ~meta_effect_sizes-1.96*meta_effect_sizes_se,
              ymax = ~meta_effect_sizes+1.96*meta_effect_sizes_se) %>%
  layout(showlegend = TRUE)

subplot(c1,e1, nrows = 1)

r <- plot_ly(check_data,
              y = ~agemid, 
              x = ~as.numeric(year), 
              z = ~rate, 
              color = ~year,
              type = 'scatter3d', mode = 'lines') 

mr <- plot_ly(check_data,
             y = ~agemid, 
             x = ~as.numeric(year), 
             z = ~meta_rate, 
             color = ~year,
             type = 'scatter3d', mode = 'lines') 

l <- plot_ly(check_data,
              y = ~agemid, 
              x = ~as.numeric(year), 
              z = ~logit, 
              color = ~year,
              type = 'scatter3d', mode = 'lines')

ml <- plot_ly(check_data,
        y = ~agemid, 
        x = ~as.numeric(year), 
        z = ~meta_logit, 
        color = ~year,
        type = 'scatter3d', mode = 'lines') 


subplot(r,mr) %>% layout(showlegend = FALSE)

subplot(l,ml) %>% layout(showlegend = FALSE)

plot_ly(check_data,
        y = ~agemid, 
        x = ~as.numeric(year), 
        z = ~meta_effect_sizes, 
        color = ~year,
        type = 'scatter3d', mode = 'lines') 

