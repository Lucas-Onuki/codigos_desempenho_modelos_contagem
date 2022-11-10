######################################################## Pacotes utilizados ----

#para os ajustes referência 
if(!require(galmss)) install.packages("gamlss")
library(gamlss)

#para o ajuste bell
if(!require(bellreg)) install.packages("bellreg")
library(bellreg)

#para o ajuste TwP
if(!require(mcglm)) install.packages("mcglm")
library(mcglm)

#para os dados da poisson tweedie ####-
if(!require(poistweedie)) install.packages("poistweedie")
library(poistweedie)

#cálculo de MSE
if(!require(MLmetrics)) install.packages("MLmetrics")
library(MLmetrics)

#para salvar o arquivo em .xlsx
if(!require(writexl)) install.packages("writexl")
library(writexl)

###################################################### Construíndo os dados ----

#tamanho da amostra
n = 500

#indicadores - sequência
indicadores = seq(from = 1, to = 100, by = 1)

#definição dos parâmetros
beta0 = 0.5
beta1 = 1
beta2 = -1
betas = c(beta0,beta1,beta)

#definição das variáveis 
x1 <- list()
x2 <- list()

for (i in indicadores) {
  
  set.seed(i)
  
  x1[[i]] = rnorm(n)
  x2[[i]] = runif(n)
}


#vetor de médias
mu <- list()

for (i in indicadores) {
  
  set.seed(i)
  
  mu[[i]] <- exp(beta0+beta1*x1[[i]]+beta2*x2[[i]])
}

library(LambertW)
theta <- list()

for (i in indicadores) {
  
  theta[[i]] = W(mu[[i]])
  
}

#variável resposta
y <- list() #iniciando nossa lista
for (i in indicadores) {
  y[[i]] <- rep(0,500) #criando uma lista com zeros, que será modificada abaixo
}

for (i in indicadores) {
  set.seed(i)
  for (j in 1:n) {
    y[[i]][j] <- rbell(1,theta[[i]][j]) #modificando os valores da lista
  }
}


#construindo o conj. de dados
dt <- list()
for (i in indicadores) {
  dt[[i]] <- as.data.frame(cbind(y[[i]],x1[[i]],x2[[i]]))
}

# legenda
#V1 - "y"
#v2 - "x1"
#v3 - "x2"

###################################################### Análise Frequentista ----

# Poisson - via gamlss: ----

p.gamlss <- NULL

#guardando ajustes em lista
for (i in 1:100) {
  p.gamlss[[i]] <- gamlss(dt[[i]]$V1 ~ dt[[i]]$V2 + dt[[i]]$V3,
                          family = PO, 
                          data = dt[[i]])
}

#salvando os coeficientes
p.coefs <- NULL
for (i in 1:100) {
  p.coefs[[i]] <- p.gamlss[[i]]$mu.coefficients
}

#manipulação do dataframe com os coeficientes
p.coefs_save <- as.data.frame(p.coefs)
p.coefs_save <- t(p.coefs_save) 
p.coefs_save <- as.data.frame(p.coefs_save)

#cálculo dos MSE's:
MSE_intercept <- MSE(p.coefs_save$`(Intercept)`,0.5)
#0.2555835
MSE_x1 <- MSE(p.coefs_save$`dt[[i]]$V2`,1)
#2.088576
MSE_x2 <- MSE(p.coefs_save$`dt[[i]]$V3`,-1)
#0.01247919

#salvando os valores dos coeficientes

write_xlsx(p.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_poisson_Bell.xlsx")

# Binomial Negativa - via gamlss: ----

#guardando ajustes em lista

nb.gamlss <- NULL

for (i in 1:100) {
  nb.gamlss[[i]] <- gamlss(dt[[i]]$V1 ~ dt[[i]]$V2 + dt[[i]]$V3,
                           family = NBF, 
                           data = dt[[i]])
}

nb.gamlss[[1]]$mu.coefficients
theta = 1/exp(nb.gamlss[[1]]$sigma.coefficients)

#salvando os coeficientes
nb.coefs <- NULL
nb.theta <- NULL
for (i in 1:100) {
  nb.coefs[[i]] <- nb.gamlss[[i]]$mu.coefficients
  nb.theta[[i]] <- 1/exp(nb.gamlss[[i]]$sigma.coefficients)
}

#manipulação do dataframe com os coeficientes
nb.coefs_save_1 <- as.data.frame(nb.coefs)
nb.coefs_save_2 <- as.data.frame(nb.theta)
nb.coefs_save_1 <- t(nb.coefs_save_1) 
nb.coefs_save_2 <- t(nb.coefs_save_2) 
nb.coefs_save_1 <- as.data.frame(nb.coefs_save_1)
nb.coefs_save_2 <- as.data.frame(nb.coefs_save_2)

nb.coefs_save <- cbind(nb.coefs_save_1,nb.coefs_save_2)

names(nb.coefs_save)[4] <- "theta"

#cálculo dos MSE's:
MSE_intercept_nb <- MSE(nb.coefs_save$`(Intercept)`,0.5)
#0.00437283
MSE_x1_nb <- MSE(nb.coefs_save$`dt[[i]]$V2`,1)
#0.008870212
MSE_x2_nb <- MSE(nb.coefs_save$`dt[[i]]$V3`,-1)
#0.0123525
MSE_theta_nb <- MSE(nb.coefs_save$theta,2)
#0.4159196

#salvando os valores dos coeficientes

write_xlsx(nb.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_NB_Bell.xlsx")




# Bell - via bellreg: ----

#guardando ajustes em lista

bell.bellreg <- NULL

for (i in 1:100) {
  bell.bellreg[[i]] <- bellreg::bellreg(dt[[i]]$V1 ~ dt[[i]]$V2 + dt[[i]]$V3,
                                        data = dt[[i]],
                                        approach = "mle")
}

coef(bell.bellreg[[1]])

#salvando os coeficientes
bell.coefs <- NULL

for (i in 1:100) {
  bell.coefs[[i]] <- coef(bell.bellreg[[i]])
}

#manipulação do dataframe com os coeficientes
bell.coefs_save <- as.data.frame(bell.coefs)
bell.coefs_save <- t(bell.coefs_save) 
bell.coefs_save <- as.data.frame(bell.coefs_save)

#cálculo dos MSE's:
MSE_intercept_b <- MSE(bell.coefs_save$`(Intercept)`,0.5)
#0.004398341
MSE_x1_b <- MSE(bell.coefs_save$`dt[[i]]$V2`,1)
#0.008808564
MSE_x2_b <- MSE(bell.coefs_save$`dt[[i]]$V3`,-1)
#0.01241896

#salvando os valores dos coeficientes

write_xlsx(bell.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_bell_Bell.xlsx")


# Tweedie Poisson - via mcglm: ----

#guardando ajustes em lista

twp.mcglm <- NULL
for (i in 1:100) {
  twp.mcglm[[i]] <- mcglm::mcglm(linear_pred = c(dt[[i]]$V1 ~ dt[[i]]$V2 + dt[[i]]$V3),
                                 matrix_pred = list(mc_id(data = dt[[i]])),
                                 link = "log",
                                 variance = "poisson_tweedie",
                                 power_fixed = FALSE,
                                 data = dt[[i]],
                                 control_algorithm = list(verbose = FALSE,
                                                          max_iter = 100,
                                                          tunning = 0.5,
                                                          correct = FALSE))
}

twp.mcglm[[1]]$Regression
twp.mcglm[[1]]$Covariance[2] #parâmetro de dispersão - lambda

#salvando os coeficientes
twp.coefs <- NULL
twp.lambda <- NULL
for (i in 1:100) {
  twp.coefs[[i]] <- twp.mcglm[[i]]$Regression
  twp.lambda[[i]] <- twp.mcglm[[i]]$Covariance[2]
}

#manipulação do dataframe com os coeficientes
twp.coefs_save_1 <- as.data.frame(twp.coefs)
twp.coefs_save_2 <- as.data.frame(twp.lambda)
twp.coefs_save_1 <- t(twp.coefs_save_1) 
twp.coefs_save_2 <- t(twp.coefs_save_2) 
twp.coefs_save_1 <- as.data.frame(twp.coefs_save_1)
twp.coefs_save_2 <- as.data.frame(twp.coefs_save_2)

twp.coefs_save <- cbind(twp.coefs_save_1,twp.coefs_save_2)

names(twp.coefs_save)[4] <- "lambda"

#cálculo dos MSE's:
(MSE_intercept_twp <- MSE(twp.coefs_save$V1,0.5))
#0.004388261
(MSE_x1_twp <- MSE(twp.coefs_save$V2,1))
#0.008866746
(MSE_x2_twp <- MSE(twp.coefs_save$V3,-1))
#0.01260233
(MSE_lambda_twp <- MSE(twp.coefs_save$lambda,2)) #???
#0.4159196

#salvando os valores dos coeficientes

write_xlsx(twp.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_twp_Bell.xlsx")


###############################################################################- 
#ref.: ----

twp.mcglm <- mcglm::mcglm(linear_pred = c(y ~ x1 + x2),
                          matrix_pred = list(mc_id(data = dt)),
                          link = "log",
                          variance = "poisson_tweedie",
                          power_fixed = FALSE,
                          data = dt,
                          control_algorithm = list(verbose = FALSE,
                                                   max_iter = 100,
                                                   tunning = 0.5,
                                                   correct = FALSE))


################################################################ Referência TwP-

#m2 <- mcglm(linear_pred = c(ncap ~ est * (des + I(des^2))),
#            matrix_pred = list(mc_id(data = capdesfo)),
#            link = "log",
#            variance = "poisson_tweedie",
#            power_fixed = FALSE,
#            data = capdesfo,
#            control_algorithm = list(verbose = FALSE,
#                                     max_iter = 100,
#                                     tunning = 0.5,
#                                     correct = FALSE))

