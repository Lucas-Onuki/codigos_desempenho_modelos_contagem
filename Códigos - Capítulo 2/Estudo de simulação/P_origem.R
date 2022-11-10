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
  
  x1[[i]] = runif(n,-1,1)
  x2[[i]] = rbinom(n,1,0.5)
}


#vetor de médias
mu <- list()

for (i in indicadores) {
  
  set.seed(i)
  
    mu[[i]] <- exp(beta0+beta1*x1[[i]]+beta2*x2[[i]])
}

View(as.data.frame(mu[[1]]))

mu[[1]][1]

#variável resposta
y <- list() #iniciando nossa lista
for (i in indicadores) {
  y[[i]] <- rep(0,500) #criando uma lista com zeros, que será modificada abaixo
}


for (i in indicadores) {
  set.seed(i)
  for (j in 1:n) {
    y[[i]][j] <- rpois(1,mu[[i]][j]) #modificando os valores da lista
  }
}
View(as.data.frame(y[[1]]))

y[[1]]
#descritiva
#hist(y)
#mean(y)
#var(y)

#construindo o conj. de dados

dados <- list()
for (i in indicadores) {
  dados[[i]] <- as.data.frame(cbind(y[[i]],x1[[i]],x2[[i]]))
}

View(as.data.frame(dados[[1]]))

# legenda
#V1 - "y"
#v2 - "x1"
#v3 - "x2"

###################################################### Análise Frequentista ----

# Poisson - via gamlss: ----

p.gamlss <- NULL

#guardando ajustes em lista
for (i in 1:100) {
  p.gamlss[[i]] <- gamlss(dados[[i]]$V1 ~ dados[[i]]$V2 + dados[[i]]$V3,
                        family = PO, 
                        data = dados[[i]])
}

#salvando os coeficientes
p.coefs <- NULL
for (i in 1:100) {
  p.coefs[[i]] <- p.gamlss[[i]]$mu.coefficients
}
p.coefs[[1]]
#manipulação do dataframe com os coeficientes
p.coefs_save <- as.data.frame(p.coefs)
p.coefs_save <- t(p.coefs_save) 
p.coefs_save <- as.data.frame(p.coefs_save)

#cálculo dos MSE's:
MSE_intercept <- MSE(p.coefs_save$`(Intercept)`,0.5)
#0.2555835
MSE_x1 <- MSE(p.coefs_save$`dados[[i]]$V2`,1)
#2.088576
MSE_x2 <- MSE(p.coefs_save$`dados[[i]]$V3`,-1)
#0.01247919

#salvando os valores dos coeficientes

write_xlsx(p.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_poisson.xlsx")

# Binomial Negativa - via gamlss: ----

#guardando ajustes em lista

nb.gamlss <- list()

for (i in 1:100) {
  nb.gamlss[[i]] <- gamlss(dados[[i]]$V1 ~ dados[[i]]$V2 + dados[[i]]$V3,
                          family = NBF, 
                          data = dados[[i]])
}

dados[[5]]
nb.gamlss <- gamlss(dados[[4]]$V1 ~ dados[[4]]$V2 + dados[[4]]$V3,
                    family = NBF, 
                    data = dados[[4]])

nb.gamlss$mu.coefficients
# Bell - via bellreg: ----

#guardando ajustes em lista

bell.bellreg <- NULL

for (i in 1:100) {
  bell.bellreg[[i]] <- bellreg::bellreg(dados[[i]]$V1 ~ dados[[i]]$V2 + dados[[i]]$V3,
                                        data = dados[[i]],
                                        approach = "mle")
}

coef(bell.bellreg[[1]])

#salvando os coeficientes
bell.coefs <- NULL

for (i in 1:100) {
  bell.coefs[[i]] <- coef(bell.bellreg[[i]])
}



bell.bellreg <- bellreg::bellreg(y ~ x1 + x2, data = dados, approach = "mle")

#manipulação do dataframe com os coeficientes
bell.coefs_save <- as.data.frame(bell.coefs)
bell.coefs_save <- t(bell.coefs_save) 
bell.coefs_save <- as.data.frame(bell.coefs_save)

#cálculo dos MSE's:
MSE_intercept_b <- MSE(bell.coefs_save$`(Intercept)`,0.5)
#0.3201318
MSE_x1_b <- MSE(bell.coefs_save$`dados[[i]]$V2`,1)
#2.533345
MSE_x2_b <- MSE(bell.coefs_save$`dados[[i]]$V3`,-1)
#0.003478542

#salvando os valores dos coeficientes

write_xlsx(bell.coefs_save,"D:/Lucas Akio/.Mestrado/.Simulação/estudo_comparativo/coefs_bell.xlsx")


# Tweedie Poisson - via mcglm: ----

#guardando ajustes em lista

for (i in 1:100) {
  twp.mcglm[[i]] <- mcglm::mcglm(linear_pred = c(dados[[i]]$V1 ~ dados[[i]]$V2 + dados[[i]]$V3),
                                 matrix_pred = list(mc_id(data = dados[[i]])),
                                 link = "log",
                                 variance = "poisson_tweedie",
                                 power_fixed = FALSE,
                                 data = dados[[i]],
                                 control_algorithm = list(verbose = FALSE,
                                                          max_iter = 100,
                                                          tunning = 0.5,
                                                          correct = FALSE))
}

mcglm::mcglm(linear_pred = c(dados[[1]]$V1 ~ dados[[1]]$V2 + dados[[1]]$V3),
             matrix_pred = list(mc_id(data = dados[[1]])),
             link = "log",
             variance = "poisson_tweedie",
             power_fixed = FALSE,
             data = dados[[1]],
             control_algorithm = list(verbose = FALSE,
                                      max_iter = 100,
                                      tunning = 0.5,
                                      correct = FALSE))


#ref.:


teste <- mcglm::mcglm(linear_pred = c(data$y ~ data$x1 + data$x2),
                          matrix_pred = list(mc_id(data = data)),
                          link = "log",
                          variance = "poisson_tweedie",
                          power_fixed = FALSE,
                          data = data,
                          control_algorithm = list(verbose = FALSE,
                                                   max_iter = 10,
                                                   tunning = 0.5,
                                                   correct = FALSE))

aux <- list(mc_id(data = data))

aux[[1]]$Z0
library(matrixcalc)
is.positive.definite(as.matrix(aux[[1]]$Z0), tol=1e-8)

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