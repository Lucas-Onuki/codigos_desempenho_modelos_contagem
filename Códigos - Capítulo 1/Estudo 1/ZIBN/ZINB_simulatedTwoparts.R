#################################################################################################### pacotes necessários

#para contar o tempo de execu??o
library(tictoc) 

#an?lise cl?ssica
library(glmmADMB)
library(gamlss)

#an?lise bayesiana
library(INLA) 
#library(MCMCglmm)
#library(arm)
library(runjags)
library(rjags)
library(R2jags)
library(digest)
library(R2admb)
library(ggplot2)
library(pscl)
library(boot)
library(MASS)
#install.packages("brms")
library(brms)

####################################################################### Fixando o diretório e gerando os dados simulados

set.seed(1234) #fixando a semente

tic("Gerando dos dados")

#tamanho da amostra
n <- 500

#definição dos parâmetros
beta0 = 0.5
beta1 = 1
beta = -1



#definição das variáveis 
x1 <- rnorm(n,0,1)
x2 <- rnorm(n,0,1)
#x2 <- runif(n)
mu = exp(beta0 + beta1*x1 + beta*x2)
summary(mu)
#size
phi=2

#####################################

#probability of zeros
delta0 = -4
delta1 = 0.5
delta2 = -3.5
pi<-plogis(delta0+delta1*x1+delta2*x2)
summary(pi)

?plogis()

y<-NULL
for (i in 1:n) {
  y.zero <- rbinom(n, prob =  1-pi[i], size = 1)
  y[i] <- ifelse(y.zero == 0, 0, rnbinom(n = n, mu = mu[i], size = phi))
}
table(y.zero)

#descritiva
hist(y)
mean(y)
var(y)
table(y)


#construindo o conj. de dados
data<-cbind(y,x1,x2)
data<-as.data.frame(data)
#View(data)
toc()

##INFERENCIA CLASSICA
######################################################################################## ANÁLISE FREQUENTISTA (cLÁSSICA)
#classico
#pscl
if(!require(pscl)) install.packages("pscl")  #ok
library(pscl)

tic("Ajuste utilizando o pacote pscl")
zinb.zeroinfl <- zeroinfl(y ~ x1 + x2 | x1 + x2,  data = data, dist = "negbin")
summary(zinb.zeroinfl)
toc()
AIC(zinb.zeroinfl)
#size estimado
zinb.zeroinfl$theta

########################################################################################################################
#classico
#gamlss
if(!require(gamlss)) install.packages("gamlss")  #ok
library(gamlss)

tic("Ajuste utilizando pacote gamlss")
zinb.gamlss<-gamlss(formula=y ~ x1 +x2, nu.formula =y ~ x1 +x2, family =ZINBI, data =data,
                                      trace =TRUE)
summary(zinb.gamlss)
toc()
#size estimado
sigma<-exp(zinb.gamlss$sigma.coefficients)
1/sigma



########################################################################################################################
#VGAM
#classico
if(!require(VGAM)) install.packages("VGAM");library(VGAM)

tic("Ajuste utilizando pacote VGAM") 
zinb.VGAM <- vglm(y ~ x1 +x2, zinegbinomial(zero=1), data=data)
summary(zinb.VGAM)
toc()
AIC(zinb.VGAM)
zinb.VGAM
#Nao trabalha bem

########################################################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB)

tic("Ajuste utilizando pacote glmmTMB") 
zinb.glmmTMB <- glmmTMB(y ~ x1 +x2, ziformula=~x1+x2, family=nbinom2, data=data)
summary(zinb.glmmTMB)
toc() 
#Overdispersion parameter for nbinom2 family (): 2.74

## INFERENCIA BAYESIANA
###################################################################################################### ANÁLISE BAYESIANA
#bayesiano
#instalar previamente o JAGS
#JAGS 
if(!require(runjags)) install.packages("runjags")
library(runjags)

library("R2jags")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tic("Ajuste utilizando o pacote jags")
datas <- list ("y"= y,"x1"=x1,"x2"=x2,"n"=n)  
inits <- function(){list(beta=c(0.5,1,-1), alpha = c(-4.5,0.5,-4.5), phi=2)} 
param<-c("beta","alpha","phi")
#zinb.jags<-run.jags(model="modelzinb.txt",monitor= param,data= datas,inits= inits,n.chains=3,
#                   sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)
#summary(zinb.jags)
#extract(zinb.jags, what='dic')
zinb.jags<-jags(model.file="modelzinb.txt",parameters.to.save= param,data= datas,inits= inits,n.chains=3,
                   n.iter=10000,n.burnin=1000)
zinb.jags.upd <- update(zinb.jags, n.iter = 1000)
print(zinb.jags)

toc()
plot(zinb.jags)
traceplot(zinb.jags)



########################################################################################################################
#bayesiano
#MCMCglmm

if(!require(MCMCglmm)) install.packages("MCMCglmm")
#MCMCglmm nao trabalha ZINB



########################################################################################################################
#bayesiano
#instalar previamente stan
#brms
if(!require(stats)) install.packages("stats")
library(stats)


if(!require(brms)) install.packages("brms")
library(brms)

tic("Ajuste utilizando pacote brms")
zinb.brm <- brm(bf(y ~ x1 + x2, zi ~ x1+x2), 
               data = data, family = zero_inflated_negbinomial(),
               prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)
summary(zinb.brm)
toc()
WAIC(zinb.brm)


#bayesiano
#glmmADMB
#nao trabalha 2 partes

########################################################################################################################
#bayesiano
#inla
#nao trabalha 2 partes
########################################################################################################################






