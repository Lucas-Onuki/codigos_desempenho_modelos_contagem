#################################################################################################### pacotes necess?rios

#para contar o tempo de execução
library(tictoc) 

#análise clássica
library(glmmADMB)

install.packages("R2admb", lib="D:/Repositorio - R")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source", lib="D:/Repositorio - R")

library(gamlss)

#análise bayesiana
library(INLA)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#library(MCMCglmm)
#library(arm)
#install.packages("~/Downloads/brms_2.16.3.tar")
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

############################### Fixando o diretório e gerando os dados simulados

set.seed(1) #fixando a semente

tic("Gerando dos dados")

#tamanho da amostra
n = 500

#definição dos parâmetros
beta0 = 0.5
beta1 = 1
beta = -1

#definição das variáveis 
x1 = runif(n,-1,1)
x2 = rbinom(n,1,0.5)
size = 2

#vetor de médias
mu <- exp(beta0 + beta1*x1 + beta*x2)
mu

#variável resposta
y<-NULL
for (i in 1:n) {
  y[i]<-rnbinom(1,mu = mu[i], size = size)  
}
y

#descritiva
hist(y)
mean(y)
var(y)

#construindo o conj. de dados
data<-cbind(y,x1,x2)
data<-as.data.frame(data)
toc()

######################################################################################## AN?LISE FREQUENTISTA (cL?SSICA)

tic("Ajuste utilizando o pacote base")

#base
nb.glm <- glm.nb(y ~ x1 +x2, data = data)
summary(nb.glm)

toc()

########################################################################################################################
#clasico
#mgcv
if(!require(mgcv)) install.packages("mgcv")
library(mgcv)
tic("Ajuste utilizando pacote mgcv")
nb.mgcv = gam(y~x1+x2,family=nb(),data=data)
summary(nb.mgcv)
AIC(nb.mgcv)
toc()

########################################################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB)

tic("Ajuste utilizando pacote glmmTMB") 
nb.glmmTMB<-glmmTMB(y ~ x1 +x2, family=nbinom1(link = "log"), data=data,)
summary(nb.glmmTMB)
toc() 

1/exp(0.582)

########################################################################################################################

tic("Ajuste utilizando pacote glmmADMB")
library(glmmADMB)
#glmmADMB
nb.glmmADMB <- glmmadmb(y ~ x1 +x2, family="nbinom",data=data)
summary(nb.glmmADMB)
#Observa??o: No help tem um exemplo para simular dados zip poisson
nb.glmmADMB$phi
toc()

########################################################################################################################

tic("Ajuste utilizando pacote gamlss")
library(gamlss)
#gamlss
nb.gamlss <- gamlss(y ~ x1 +x2, family=NBF, data=data) 
summary(nb.gamlss)
nb.gamlss$sigma.coefficients
1/exp(-0.9614)
toc()

###################################################################################################### AN?LISE BAYESIANA

### arm n?o funcionou ###

tic("Ajuste utilizando pacote brms")

#brms
nb.brm <- brm(y ~ x1 + x2, family = "negbinomial", data = data)
summary(nb.brm)
WAIC(nb.brm)

toc()

########################################################################################################################

tic("Ajuste utilizando o pacote jags")

#JAGS
datas <- list ("y"=y,"x1"=x1,"x2"=x2,"n"=n)  
inits <- function(){list(beta=c(0.5,1,-1))}
param <- c("beta","deviance")
nb.jags <- run.jags(model="modelnb.txt",monitor= param,data= datas,inits= inits,n.chains=3,
                 sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)

#print(p.jags)
summary(nb.jags)
#summary(p.jags)[,4]
extract(nb.jags, what='dic')

toc()

########################################################################################################################

tic("Ajuste utilizando pacote inla")

#inla
nb.inla <- inla(y~x1+x2, family='nbinomial',
               data=data,control.family=list(link='log'),
               control.predictor=list(link=1, compute=TRUE),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(nb.inla)

toc()

########################################################################################################################

tic("Ajuste utilizando pacote brms")

#brms
nb.brm <- brm(y ~ x1 + x2, family = "negbinomial", data = data)
summary(nb.brm)
WAIC(nb.brm)

toc()

######################################################### Cálculo do afastamento

#CADA MODELO E DIFERENTE, DAR UMA OLHADA NA CHAMADA DOS COEFICIENTES

vetor_de_afastamento <- NULL

#FREQUENTISTA (cLÁSSICA)
(afastamento_glm <- (((0.5 - nb.glm$coefficients[1])^2 + (1 - nb.glm$coefficients[2])^2 + ((-1) - nb.glm$coefficients[3])^2)/3))
vetor_de_afastamento[1] <- afastamento_glm
(afastamento_glmmADMB <- (((0.5 - nb.glmmADMB$b[1])^2 + (1 - nb.glmmADMB$b[2])^2 + ((-1) - nb.glmmADMB$b[3])^2)/3))
vetor_de_afastamento[2] <- afastamento_glmmADMB
(afastamento_gamlss <- (((0.5 - nb.gamlss$mu.coefficients[1])^2 + (1 - nb.gamlss$mu.coefficients[2])^2 + ((-1) - nb.gamlss$mu.coefficients[3])^2)/3))
vetor_de_afastamento[3] <- afastamento_gamlss
  
#ANÁLISE BAYESIANA                                                                                                                                
(afastamento_jags <- (((0.5 - 0.4555522)^2 + (1 - 0.8338336)^2 + ((-1) - (-0.9656084))^2)/3)) 
vetor_de_afastamento[4] <- afastamento_jags
(afastamento_inla <- (((0.5 - 0.45800)^2 + (1 - 0.83200)^2 + ((-1) - (-0.96500))^2)/3))  
vetor_de_afastamento[5] <- afastamento_inla
(afastamento_brm <- (((0.5 - 0.40)^2 + (1 - 1.01)^2 + ((-1) - (-0.95))^2)/3)) 
vetor_de_afastamento[6] <- afastamento_brm

vetor_de_afastamento
