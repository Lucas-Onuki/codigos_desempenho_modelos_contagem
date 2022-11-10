#################################################################################################### pacotes necess?rios

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

####################################################################### Fixando o diret?rio e gerando os dados simulados

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

#size
phi=2

#probability of zeros
pi<-0.2

y<-NULL
for (i in 1:n) {
  y.zero<-rbinom(n, prob =  1-pi, size = 1)
  y[i]<-ifelse(y.zero == 0, 0, rnbinom(n = n, mu = mu[i], size = phi))
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
#gamlss
if(!require(gamlss)) install.packages("gamlss") #ok
library(gamlss)

tic("Ajuste utilizando pacote gamlss")
zinb.gamlss<-gamlss(formula=y ~ x1 +x2, family =ZINBI, data =data, method =RS(100),
                    trace =TRUE)
summary(zinb.gamlss)
toc()

#size estimado
sigma<-exp(zinb.gamlss$sigma.coefficients)
1/sigma


#probabilidade de zero estimada
nu<-zinb.gamlss$nu.coefficients
exp(nu)/(1+exp(nu))


########################################################################################################################

#https://uvastatlab.github.io/2019/08/29/simulating-data-for-count-models/
#VGAM
#classico
if(!require(VGAM)) install.packages("VGAM");library(VGAM)                                           #erro na função
#(zero = 1)
tic("Ajuste utilizando pacote VGAM") 
zinb.VGAM <- vglm(y ~ x1 +x2, zinegbinomialff(zero = c(2,3)) 
                 ,data=data)
summary(zinb.VGAM)
toc()
AIC(zinb.VGAM)
#incorreto

########################################################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB) #ok na saída da binomial usa logit como no gamlss

tic("Ajuste utilizando pacote glmmTMB") 
zinb.glmmTMB<-glmmTMB(y ~ x1 +x2, ziformula=~1, family=nbinom2, data=data)
summary(zinb.glmmTMB)
toc() 
#size estimado
exp(zinb.glmmTMB$fit$par[5])
#probabilidade de zero estimada
plogis(zinb.glmmTMB$fit$par[4])

#clasico
#mgcv
#mgcv nao trabalha ZINB


## INFERENCIA BAYESIANA
###################################################################################################### ANÁLISE BAYESIANA
#bayesiano
#instalar previamente o JAGS
#JAGS 
if(!require(runjags)) install.packages("runjags")                                               #ver
library(runjags)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tic("Ajuste utilizando o pacote jags")
datas <- list ("y"= y,"x1"=x1,"x2"=x2,"n"=n)  
#inits <- function(){list(beta=c(0.5,1,-1), q = 0.2, phi=2)}
inits <- function(){list(beta=c(0.5,1,-1), q = 0.2)}
param<-c("beta","q","phi")
zinb.jags<-run.jags(model="modelzinb1.txt",monitor= param,data= datas,inits= inits,n.chains=3,
                   sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)
summary(zinb.jags)
toc()
extract(zinb.jags, what='dic')

#MCMCglmm nao trabalha ZINB


########################################################################################################################
#bayesiano
#instalar previamente stan
#brms
if(!require(stats)) install.packages("stats")                                 #ok
library(stats)

tic("Ajuste utilizando pacote brms")                                                 
zinb.brm <- brm(bf(y ~ x1 + x2), 
               data = data, family = zero_inflated_negbinomial(),
prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)

summary(zinb.brm)
WAIC(zinb.brm)
toc()


######################################################################################################
#bayesiano
#glmmADMB
if(!require(R2admb)) install.packages("R2admb");library(R2admb)                 #ok
if(!require(glmmADMB)) 
  install.packages("glmmADMB", 
                   repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                           getOption("repos")),
                   type="source")
library(glmmADMB)

tic("Ajuste utilizando pacote glmmADMB")
zinb.glmmADMB<-glmmadmb(y ~ x1 + x2, family="nbinom",data=data, zeroInflation = TRUE)
summary(zinb.glmmADMB)
toc()
#size estimado
zinb.glmmADMB$alpha
#probabilidade de zero estimada
zinb.glmmADMB$pz

########################################################################################################################
#bayesiano
#inla
if(!require(INLA))                                                                      #ok
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla -download.org/R/stable"), dep=TRUE) 

library(INLA) 

tic("Ajuste utilizando pacote inla")
zinb.inla = inla(y~x1+x2, family = "zeroinflatednbinomial0", 
                data = data,control.family=list(link='log'), 
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE)) 
summary(zinb.inla)
toc()


########################################################################################################################


