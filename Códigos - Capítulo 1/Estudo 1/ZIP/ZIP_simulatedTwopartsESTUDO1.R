#################################################################################################### pacotes necessários
#ZIP two parts
###############
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
alpha <- c(0.4, -0.5, 0.9)

x1 <- rnorm(n,0,1)
x2 <- rnorm(n,0,1)
mu <- exp(beta[1]+beta[2]*x1+beta[3]*x2)
pi <- plogis(alpha[1]+alpha[2]*x1+alpha[3]*x2)

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
#pscl
if(!require(pscl)) install.packages("pscl")
library(pscl)

tic("Ajuste utilizando o pacote pscl")
zip.zeroinfl <- zeroinfl(y ~ x1 + x2,  data = data, dist = "poisson", EM = T)
summary(zip.zeroinfl)
AIC(zip.zeroinfl)
toc()


########################################################################################################################
#classico
#gamlss
if(!require(gamlss)) install.packages("gamlss")
library(gamlss)

tic("Ajuste utilizando pacote gamlss")
zip.gamlss<-gamlss(formula=y ~ x1 +x2, sigma.formula =y ~ x1 +x2, family =ZIP, data =data, method =RS(100),
                   trace =TRUE)
summary(zip.gamlss)
toc()

########################################################################################################################
#VGAM
#classico
if(!require(VGAM)) install.packages("VGAM");library(VGAM)

tic("Ajuste utilizando pacote VGAM") 
zip.VGAM <- vglm(y ~ x1 +x2, zipoisson, data=data)
summary(zip.VGAM)
toc()
AIC(zip.VGAM)

########################################################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB)

tic("Ajuste utilizando pacote glmmTMB") 
zip.glmmTMB <- glmmTMB(y ~ x1 +x2, ziformula=~x1+x2, family=poisson, data=data)
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
inits <- function(){list(beta=c(0.5,1,-1), alpha = c(0.4,-0.5,0.9),g=rep(1,n))}
param<-c("beta","deviance","alpha")
zip.jags<-run.jags(model="modelzip.txt",monitor= param,data= datas,inits= inits,n.chains=3,
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
zip.MCMCglmm <- MCMCglmm(y~trait+trait:x1+trait:x2,  rcov = ~trait:units, data = data,
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
zip.brm <- brm(bf(y ~ x1 + x2, zi ~ x1+x2), 
                   data = data, family = zero_inflated_poisson())
summary(zip.brm)
WAIC(zip.brm)
toc()


