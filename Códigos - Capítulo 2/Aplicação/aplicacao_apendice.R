####################################################### pacotes necessários ----

require(pscl)                                          # para obtenção dos dados
require(gamlss)                                  # para o ajuste do modelo ZIBN2
library(bellreg)                                             #para o ajuste Bell
library(mcglm)                                                #para o ajuste TwP

###############################################################################-

########################################################## leitura de dados ----

# limpando o console
rm(list = ls(all.names = TRUE)) 

# leitura do banco
data("bioChemists")
dados = attach(bioChemists)
dados = as.data.frame(bioChemists)
# variável resposta - NÚMERO DE ARTIGOS PRODUZÍDOS NOS 3 ÚLTIMOS ANOS DE PHD
hist(art,
     ylab="Histogram",
     col = "green",
     main="Histograma para variável resposta")

###############################################################################-

################################ MODELAGEM #####################################

################################ modelo Zero-inflacionado Binomial Negativo ----

# pacote utilizado gamlss
zinb.gamlss2 <- gamlss::gamlss(formula = art ~ fem + kid5 + ment,
                               nu.formula = art ~ ment,
                               family = ZINBI,
                               data = bioChemists,
                               trace = TRUE)

summary(zinb.gamlss2) # sumário do modelo

# size estimado
sigma <- exp(zinb.gamlss2$sigma.coefficients)
1/sigma

#AIC
zinb.gamlss2$aic #computando o AIC

###############################################################################-

############################################################### modelo Bell ----

# pacote utilizado bellreg
bell.zibellreg <- bellreg::zibellreg(formula = ment |art ~ fem + kid5 + ment,
                          data = dados,
                          approach = "mle",
                          hessian = T)



summary(bell.zibellreg) # sumário do modelo

#AIC
bell.zibellreg$AIC #computando o AIC

###############################################################################-

#################################################### modelo Tweedie-Poisson ----

twp.mcglm <- mcglm::mcglm(linear_pred = c(art ~ fem + kid5 + ment),
                          matrix_pred = list(mc_id(data = bioChemists)),
                          link = "log",
                          variance = "poisson_tweedie",
                          power_fixed = FALSE,
                          data = bioChemists,
                          control_algorithm = list(verbose = FALSE,
                                                   max_iter = 100,
                                                   tunning = 0.5,
                                                   correct = FALSE))

summary(twp.mcglm) # sumário do modelo

gof(twp.mcglm) # o critério de informação aqui é o pAIC p de pseudo, essa medida
# pelo que li é equivalente ao AIC normal, porém utilizada para
# modelos com pseudo-logverossimilhança

###############################################################################-