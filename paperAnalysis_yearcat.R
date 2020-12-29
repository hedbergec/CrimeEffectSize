rm(list = ls(all.names = TRUE))

library(data.table)
library(plotly)
library(glmnet)
library(reticulate)

sklearn <- import("sklearn")
np <- import("numpy")
rand <- import("random")

load("dashData.Rdata")

set.seed(101)

yvar <- "othassault"
yvar <- "totaldrug"
yvar <- "aggassault"

poly <- 2

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
  filter(is.infinite(!!sym(paste0("l_",yvar,"_tot"))) == FALSE) %>%
  select(data_vars) %>%
  rename(lnagemid1 = lnagemid) %>%
  mutate(lnyear = log(as.numeric(year))-log(2015), lnpop = log(pop_tot))

list_of_states <- rownames(table(input_data$state))
list_of_years <- rownames(table(input_data$year))

lnagevals <- as.numeric(rownames(table(input_data$lnagemid1)))
agevals <- as.numeric(rownames(table(input_data$agemid)))

x <- cbind(input_data[,"lnagemid1"])

#### effects of age ####

for (p in 2:poly) {
  x[,paste0("lnagemid",p)] <- input_data[,"lnagemid1"]^p
}

#### effects of state ####

for (s in 2:length(list_of_states)) {
  x[,paste0(list_of_states[s])] <- 0
  x[which(input_data[,"state"]==list_of_states[s]),list_of_states[s]] <- 1
  x[which(input_data[,"state"]==list_of_states[1]),list_of_states[s]] <- -1
}


#### effects of year ####

for (y in 2:length(list_of_years)) {
  x[,paste0(list_of_years[y])] <- 0
  x[which(input_data[,"year"]==list_of_years[y]),list_of_years[y]] <- 1
  x[which(input_data[,"year"]==list_of_years[1]),list_of_years[y]] <- -1
}

#### center ####

for (v in c(list_of_states[-1],list_of_years[-1])){
  #x[,v] <- x[,v] - mean(x[,v])
}

#### outcome ####

out <- cbind(input_data[,paste0("l_",yvar,"_tot")])

##### Run OLS Model ####

OLS_data <- cbind(out,x)

colnames(OLS_data) <- c("y",paste0("x",seq(1:(ncol(x)))))

form <- "y~1"

for (i in seq(1:(ncol(x)))) {
  form <- paste0(form,"+x",i)
}

ols_mod <- lm(as.formula(form),
                 data = OLS_data
)

summary(ols_mod)$r.squared

ols_coefs <- t(coef(ols_mod))

ols_hat <- rep(ols_coefs[1], length(lnagevals))
for (p in 1:poly) {
  ols_hat <- ols_hat + ols_coefs[p+1]*lnagevals^p
}

ols_effect_sizes <- rep(0,length(agevals))
for (p in 1:poly) {
  ols_effect_sizes <- ols_effect_sizes +
    p*ols_mod$coefficients[p+1]*lnagevals^(p-1)/agevals
}

ols_effect_sizes <- ols_effect_sizes*sqrt(3)/pi

##### Run LASSO Model ####

lasso_mod <- glmnet(as.matrix(x),
                    as.matrix(out))

lasso_coefs <- t(coef(lasso_mod, s=min(lasso_mod$lambda)))

lasso_hat <- rep(lasso_coefs[1], length(lnagevals))
for (p in 1:poly) {
  lasso_hat <- lasso_hat + 
    lasso_coefs[p+1]*lnagevals^p
}

lasso_effect_sizes <- rep(0,length(agevals))
for (p in 1:poly) {
  lasso_effect_sizes <- lasso_effect_sizes +
    p*lasso_coefs[p+1]*lnagevals^(p-1)/agevals
}

lasso_effect_sizes <- lasso_effect_sizes*sqrt(3)/pi

##### Run RANSAC Model (from python) ####

ransac <- sklearn$linear_model$RANSACRegressor()

ransac$random_state <- np$int(101)

ransac_mod <- ransac$fit(X = x, y = out)

ransac_coefs <- c(ransac_mod$estimator_$intercept_,ransac_mod$estimator_$coef_)

ransac_hat <- rep(ransac_coefs[1], length(lnagevals))
for (p in 1:poly) {
  ransac_hat <- ransac_hat + 
    ransac_coefs[p+1]*lnagevals^p
}

ransac_effect_sizes <- rep(0,length(agevals))
for (p in 1:poly) {
  ransac_effect_sizes <- ransac_effect_sizes +
    p*ransac_coefs[p+1]*lnagevals^(p-1)/agevals
}

ransac_effect_sizes <- ransac_effect_sizes*sqrt(3)/pi


#### Run Meta Analysis on OLS models ####

b <- c()
V <- c()
for (y in list_of_years) {
  for (s in 1:length(list_of_states)) {
    model_data <- input_data %>%
      filter(state == list_of_states[s]) %>%
      filter(year == y)
    if (nrow(model_data) > poly+1) {
      model <- lm(as.formula(paste0(paste0("l_",yvar,"_tot"),
                                    "~poly(lnagemid1,poly, raw = TRUE)")), 
                  data = model_data)
      b_set <- c()
      V_set <- c()
      for (i in 0:poly) {
        b_set[paste0("b",i)]<-coef(model)[i+1]
        for (j in 0:poly) {
          if (j >= i) {
            V_set[paste0("V_",i,"_",j)] <- vcov(model)[i+1,j+1]
          }
        }
      }
      b <- rbind(b, b_set)
      V <- rbind(V, V_set)
    }
    
  }
}


meta_mod <- mvmeta(b,V)

meta_coefs <- coef(meta_mod)

meta_hat <- rep(meta_coefs[1], length(lnagevals))
for (p in 1:poly) {
  meta_hat <- meta_hat + meta_coefs[p+1]*lnagevals^p
}

meta_effect_sizes <- rep(0,length(agevals))
for (p in 1:poly) {
  meta_effect_sizes <- meta_effect_sizes +
    p*meta_coefs[p+1]*lnagevals^(p-1)/agevals
}

meta_effect_sizes <- meta_effect_sizes*sqrt(3)/pi

#### Make Check Data ####

check_data <- input_data %>%
  group_by(agemid, year, state) %>%
  summarise(rate = mean(!!sym(paste0("r_",yvar,"_tot"))),
            logit = mean(!!sym(paste0("l_",yvar,"_tot"))))


check_data <- merge(check_data, 
                    cbind(agemid = agevals,
                          ols_logit = ols_hat,
                          ols_hat = 1e+5*exp(ols_hat)/(1+exp(ols_hat)),
                          lasso_logit = lasso_hat,
                          lasso_hat = 1e+5*exp(lasso_hat)/(1+exp(lasso_hat)),
                          ransac_logit = ransac_hat,
                          ransac_hat = 1e+5*exp(ransac_hat)/(1+exp(ransac_hat)),
                          meta_logit = meta_hat,
                          meta_hat = 1e+5*exp(meta_hat)/(1+exp(meta_hat))
                    )
)

es <- data.frame(cbind(agevals,ols_effect_sizes,lasso_effect_sizes,ransac_effect_sizes,meta_effect_sizes))

p1 <- plot_ly(check_data, 
              y =~rate, 
              x =~agemid, 
              type = "box",
              name = "Rate box plot") %>% 
  add_trace(y = ~ols_hat, 
            type = "scatter" , 
            mode = 'lines',
            name = "OLS rate") %>%
  add_trace(y = ~lasso_hat, 
            type = "scatter" , 
            mode = 'lines',
            name = "LASSO rate") %>%
  add_trace(y = ~ransac_hat, 
            type = "scatter" , 
            mode = 'lines',
            name = "RANSAC rate") %>%
  add_trace(y = ~meta_hat, 
            type = "scatter" , 
            mode = 'lines',
            name = "META rate") 

p2 <- plot_ly(check_data, 
              y =~logit, 
              x =~agemid, 
              type = "box",
              name = "Logit box plot") %>% 
  add_trace(y = ~ols_logit, 
            type = "scatter" , 
            mode = 'lines',
            name = "OLS logit") %>%
  add_trace(y = ~lasso_logit, 
            type = "scatter" , 
            mode = 'lines',
            name = "LASSO logit") %>%
  add_trace(y = ~ransac_logit, 
            type = "scatter" , 
            mode = 'lines',
            name = "RANSAC logit") %>%
  add_trace(y = ~meta_logit, 
            type = "scatter" , 
            mode = 'lines',
            name = "META logit")


p3 <- plot_ly(es, 
              y =~ols_effect_sizes, 
              x =~agevals, 
              type = "scatter" , 
              mode = 'lines',
              name = "OLS effect sizes") %>% 
  add_trace(y = ~lasso_effect_sizes, 
            type = "scatter" , 
            mode = 'lines',
            name = "LASSO effect sizes") %>%
  add_trace(y = ~ransac_effect_sizes, 
            type = "scatter" , 
            mode = 'lines',
            name = "RANSAC effect sizes") %>%
  add_trace(y = ~meta_effect_sizes, 
            type = "scatter" , 
            mode = 'lines',
            name = "META effect sizes")


subplot(p1,p2,p3, nrows = 1) %>% layout(title = yvar)
