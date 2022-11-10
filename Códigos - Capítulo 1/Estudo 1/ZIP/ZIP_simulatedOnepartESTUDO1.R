#################################################################################################### pacotes necessários
#ZIP POISSON ONE PART
#para contar o tempo de execução
if(!require(tictoc)) install.packages("tictoc");library(tictoc)


####################################################################### Fixando o diretório e gerando os dados simulados
tic("Gerando dos dados")
set.seed(1234) #fixando a semente

#tamanho da amostra e número de parâmetros
n <- 500
p <- 3

#definição dos parâmetros
beta <- c(0.5, 1, -1)
alpha <- 0.4

x1 <- rnorm(n,0,1)
x2 <- rnorm(n,0,1)
mu <- exp(beta[1]+beta[2]*x1+beta[3]*x2)
pi <- plogis(alpha)

y <- rpois(n, lambda = mu) * rbinom(n, prob = 1 - pi, size = 1)
#alternative
#zero<-rbinom(n, prob =  pi, size = 1)
#y.count <- rpois(n, lambda = mu)
#y <- y.count*(1-zero)

hist(y)
mean(y)
var(y)

x1 <- as.vector(x1)
x2 <- as.vector(x2)
#construindo o conj. de dados
data<-cbind(y,x1,x2)
data<-as.data.frame(data)
toc()

##INFERENCIA CLASSICA
######################################################################################## ANÁLISE FREQUENTISTA (cLÁSSICA)
#classico
#gamlss
if(!require(gamlss)) install.packages("gamlss")
library(gamlss)

tic("Ajuste utilizando pacote gamlss")
zip.gamlss<-gamlss(formula=y ~ x1 +x2, family =ZIP, data =data, method =RS(100),
                   trace =TRUE)
summary(zip.gamlss)
toc()

########################################################################################################################
#VGAM
#classico
if(!require(VGAM)) install.packages("VGAM");library(VGAM)

tic("Ajuste utilizando pacote VGAM") 
zip.VGAM <- vglm(y ~ x1 +x2, zipoisson(zero = 1), data=data)
summary(zip.VGAM)
AIC(zip.VGAM)
toc()


########################################################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB)

tic("Ajuste utilizando pacote glmmTMB") 
zip.glmmTMB<-glmmTMB(y ~ x1 +x2, ziformula=~1, family=poisson, data=data)
summary(zip.glmmTMB)

toc() 



## INFERENCIA BAYESIANA
###################################################################################################### ANÁLISE BAYESIANA
#bayesiano
#instalar previamente o JAGS
#JAGS 
if(!require(runjags)) install.packages("runjags")
library(runjags)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tic("Ajuste utilizando o pacote jags")
datas <- list ("y"= y,"x1"=x1,"x2"=x2,"n"=n)  
inits <- function(){list(beta=c(0.5,1,-1), alpha = 0.4,g=rep(1,n))}
param<-c("beta","deviance","alpha")
zip.jags<-run.jags(model="modelzip1.txt",monitor= param,data= datas,inits= inits,n.chains=3,
                  sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)

summary(zip.jags)
extract(zip.jags, what='dic')
toc()

########################################################################################################################
#bayesiano
#MCMCglmm
if(!require(MCMCglmm)) install.packages("MCMCglmm")
library(MCMCglmm)

tic("Ajuste utilizando pacote MCMCglm") 
zip.MCMCglmm <- MCMCglmm(y~x1+x2,  rcov = ~trait:units, data = data,
                         family = "zipoisson", verbose = FALSE,
                         nitt=4000, burnin=1000, thin=5)
summary(zip.MCMCglmm)
toc()


########################################################################################################################
#bayesiano
#instalar previamente stan
#brms
if(!require(stats)) install.packages("stats")
library(stats)


if(!require(brms)) install.packages("brms")
library(brms)

tic("Ajuste utilizando pacote brms")
zip.brm <- brm(bf(y ~ x1 + x2), 
                   data = data, family = zero_inflated_poisson())
summary(zip.brm)
WAIC(zip.brm)
toc()

######################################################################################################
#bayesiano
#glmmADMB
if(!require(R2admb)) install.packages("R2admb");library(R2admb)
if(!require(glmmADMB)) 
  install.packages("glmmADMB", 
                   repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                           getOption("repos")),
                   type="source")
library(glmmADMB)

tic("Ajuste utilizando pacote glmmADMB")
zip.glmmADMB<-glmmadmb(y ~ x1 + x2, family="poisson",data=data, zeroInflation = TRUE)
summary(zip.glmmADMB)
toc()

########################################################################################################################
#bayesiano
#inla
if(!require(INLA)) 
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla -download.org/R/stable"), dep=TRUE) 

library(INLA) 

tic("Ajuste utilizando pacote inla")
zip.inla = inla(y~x1+x2, family = "zeroinflatedpoisson0", 
                data = data,control.family=list(link='log'), 
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE)) 
summary(zip.inla)
toc()

########################################################################################################################
#clasico
#mgcv
if(!require(mgcv)) install.packages("mgcv")
library(mgcv)
tic("Ajuste utilizando pacote mgcv")
zip.mgcv = gam(y~x1+x2,family=ziP(),data=data)
summary(zip.mgcv)
AIC(zip.mgcv)
toc()



