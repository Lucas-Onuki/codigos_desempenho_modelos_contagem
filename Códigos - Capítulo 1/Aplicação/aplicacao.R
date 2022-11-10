############# Leitura dos pacotes necessários ----

# análise clássica
require(pscl)
require(gamlss)

require(package)
require(MASS)
require(nortest)
require(qcc)
require(hnp)
require(aplore3)

require(DAAG)
require(reliaR)
require(coda)
require(mcmcplots)
require(foreign)
require(psych)

require(robustbase)

# visualização
require(ggplot2)
require(car)


# análise bayesiana
require(arm)
require(INLA)

############# Leitura dos dados e descritivas ----

# limpando o console
rm(list = ls(all.names = TRUE)) 

# leitura e visualização do banco
data("bioChemists")
# data = bioChemists
# data = as.data.frame(data)
# data

#edit(bioChemists)
head(bioChemists)
attach(bioChemists)

### mais informações em: ?bioChemists ###

#View(bioChemists) #  para visualizar o banco 

# sumário do banco
summary(bioChemists)

######### DEFINICAO DE COR #####################---
color <- colorRampPalette(c("blue","lightblue"))
color(10)
################################################---

# histogramas para as variáveis resposta e explicativas (sexo e indicadora de casamento)

# variável resposta - NÚMERO DE ARTIGOS PRODUZÍDOS NOS 3 ÚLTIMOS ANOS DE PHD
hist(art,ylab="Histogram",col = color(10), main="Histograma para variável resposta")

# variáveis explicativas
ggplot(bioChemists, aes(art, fill =fem)) + geom_histogram(binwidth = 1) + 
  facet_grid(fem ~., margins = TRUE, scales = "free") +
  ggtitle("Histogramas da relação entre sexo e variável resposta") # sexo

ggplot(bioChemists, aes(art, fill =mar)) + geom_histogram(binwidth = 1) + 
  facet_grid(mar ~., margins = TRUE, scales = "free") +
  ggtitle("Histogramas da relação entre o estado civil e variável resposta") # estado civil

#----------------------------------------------#

# tabela de frequências para a variável resposta
table(art)

# Há a presença de aproximadamente 30% de zeros na variável resposta, o que é um indicador de que um ajuste
# inflacionado talvez seja o mais indicado ---

# sumário da var resposta
summary(art)

# média
mean(art)

# variância
var(art)

# São diferentes e portanto o modelo Poisson já não é o mais indicado - presença de SUPERDISPERSÃO nos dados ---

#----------------------------------------------#

# tabela de frequências para as variáveis explicativas 
table(fem)
table(mar)
table(kid5)

# boxplots de sexo e estado civil
par(mfrow=c(1,2))
boxplot(split(art, fem), col = c("lightgreen","orange"), main = "Boxplot do gênero dada a variável resposta")
boxplot(split(art, mar), col = c("lightblue","red"), main = "Boxplot do estado civil dada a variável resposta")

# boxplot ajustado elimina outliers devido a assimetria da variável resposta (Hubert and Vandervieren (2004)).
par(mfrow=c(1,2))
adjbox(split(art,fem), col = c("lightgreen","orange"), main = "Boxplot ajustado do gênero dada a variável resposta")
adjbox(split(art,mar), col = c("lightblue","red"), main = "Boxplot ajustado do estado civil dada a variável resposta")

# Os boxplots não evidenciam diferenças entre gêneros e entre estados civis, isso fica claro também via as métricas a seguir ---

with(bioChemists, tapply(art, fem, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))}))

with(bioChemists, tapply(art, mar, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

# boxplot ajustado do número de filhos
adjbox(split(art,kid5), col = c("lightgreen","orange","lightblue","red"), main = "Boxplot ajustado do número de filhos dada a variável resposta")

# Não encontramos maiores diferenças na mediana entre a variável resposta e o número de filhos, porém a variabilidade diminui à medida que o estudante tem mais filhos ---

# correlação
plot(phd,art, main = "Gráfico de pontos entre o prestígio do departamento de PHD e a variável resposta")
cor(art,phd) # 0.07333861
plot(ment,art, main = "Gráfico de pontos entre número de produções do mentor e a variável resposta")
cor(art,ment) # 0.3058618

# Encontramos una baixa correlação entre phd e a variável resposta, e uma mediana correlação entre ment e art o que indica que a maior número de publicações do orientador
# maior número de publicações do estudante indicando que esta variável pode ser uma covariável no estudo ---

#----------------------------------------------#

# As covariáveis para os ajustes são: ---

# * ment
# * mar
# * fem
# * kid5

############# Ajustes dos modelos MLG ---

# poisson ----

# gamlss
p.gamlss<-gamlss(art ~ fem + mar + kid5 + ment, family=PO, data=bioChemists) 
summary(p.gamlss)
#fit.model <- p.glm
#stepAIC(p.glm)                            # não funciona #

# AIC: 3302.349

# diagnóstico
#hnp.p.glm = hnp(p.glm, print.on=TRUE, plot=FALSE, halfnormal=F)
#plot(hnp.p.glm, las=1, pch=20, cex=1, col=c(1,1,1,2), main = "Envelope para o ajuste Poisson")
#source("http://www.ime.usp.br/~giapaula/diag_pois")

# brms
p.brm <- brm(art ~ fem + mar + kid5 + ment, family = "poisson", data = bioChemists)
summary(p.brm)
WAIC(p.brm)

# WAIC: 3318.6

# zip ----

######################################## (1) ---

#gamlss
if(!require(gamlss)) install.packages("gamlss")
library(gamlss)

zip.gamlss <- gamlss(formula=art ~ fem + mar + kid5 + ment, family =ZIP, data =bioChemists, method =RS(100),
                   trace =TRUE)
summary(zip.gamlss)

# AIC: 3253.576

#brms
if(!require(stats)) install.packages("stats")
library(stats)


if(!require(brms)) install.packages("brms")
library(brms)

zip.brm <- brm(bf(art ~ fem + mar + kid5 + ment), 
               data = bioChemists, family = zero_inflated_poisson())
summary(zip.brm)
WAIC(zip.brm)

# WAIC: 3258.7

######################################## (2) ---

#gamlss
if(!require(gamlss)) install.packages("gamlss")
library(gamlss)

zip.gamlss <- gamlss(formula=art ~ fem + mar + kid5 + ment, sigma.formula =art ~ fem + mar + kid5 + ment, family =ZIP, data =bioChemists, method =RS(100),
                   trace =TRUE)
summary(zip.gamlss)

# AIC: 3229.594

#brms
if(!require(stats)) install.packages("stats")
library(stats)
if(!require(brms)) install.packages("brms")
library(brms)

tic("Ajuste utilizando pacote brms")
zip.brm <- brm(bf(art ~ fem + mar + kid5 + ment, zi ~ fem + mar + kid5 + ment), 
               data = bioChemists, family = zero_inflated_poisson())
summary(zip.brm)
WAIC(zip.brm)

# WAIC: 3235.3 

# nb ----

# gamlss
nb.gamlss <- gamlss(art ~ fem + mar + kid5 + ment, family=NBF, data=bioChemists) 
summary(nb.gamlss)

# AIC: 3136.069

# brms
nb.brm <- brm(art ~ fem + mar + kid5 + ment, family = "negbinomial", data = bioChemists)
summary(nb.brm)
WAIC(nb.brm)

# 3134.2

# zinb ----

######################################## (1) ---

# gamlss 
if(!require(gamlss)) install.packages("gamlss")
library(gamlss)

zinb.gamlss1 <- gamlss(formula=art ~ fem + mar + kid5 + ment, family =ZIP, data = bioChemists, method =RS(100),
                   trace =TRUE)
summary(zinb.gamlss1)

# AIC: 3253.576

#brms
if(!require(stats)) install.packages("stats")                                 
library(stats)

zinb.brm1 <- brm(bf(art ~ fem + mar + kid5 + ment), 
                data = bioChemists, family = zero_inflated_negbinomial(),
                prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)

summary(zinb.brm1)
WAIC(zinb.brm1)

# WAIC: 3162.2

######################################## (2) ---

#gamlss
if(!require(gamlss)) install.packages("gamlss")  
library(gamlss)

zinb.gamlss2 <- gamlss(formula=art ~ fem + mar + kid5 + ment, nu.formula =art ~ fem + mar + kid5 + ment, family =ZINBI, data =bioChemists,
                    trace =TRUE)
summary(zinb.gamlss2)
#size estimado
sigma<-exp(zinb.gamlss2$sigma.coefficients)
1/sigma


# AIC: 3121.997

# brms 
if(!require(stats)) install.packages("stats")
library(stats)
if(!require(brms)) install.packages("brms")
library(brms)

zinb.brm2 <- brm(bf(art ~ fem + mar + kid5 + ment, zi ~  fem + mar + kid5 + ment), 
                data = bioChemists, family = zero_inflated_negbinomial(),
                prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)
summary(zinb.brm2)
WAIC(zinb.brm2)

# WAIC: 3126.2

### retirando as variáveis não significativas (MODELO FINAL) ----

#gamlss
if(!require(gamlss)) install.packages("gamlss")  
library(gamlss)

zinb.gamlss2 <- gamlss(formula=art ~ fem + kid5 + ment, nu.formula =art ~ ment, family =ZINBI, data =bioChemists,
                       trace =TRUE)
summary(zinb.gamlss2)
#size estimado
sigma<-exp(zinb.gamlss2$sigma.coefficients)
1/sigma


# AIC: 3121.997

# brms 
if(!require(stats)) install.packages("stats")
library(stats)
if(!require(brms)) install.packages("brms")
library(brms)

zinb.brm2 <- brm(bf(art ~ fem + kid5 + ment, zi ~ ment), 
                 data = bioChemists, family = zero_inflated_negbinomial(),
                 prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)
summary(zinb.brm2)
WAIC(zinb.brm2)

# WAIC: 3126.2
