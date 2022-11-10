############################################################ pacotes necessários

#para contar o tempo de execução
library(tictoc) 

#análise clássica
library(glmmADMB)
library(gamlss)

#análise bayesiana
library(INLA) 
library(MCMCglmm)
library(arm)
library(runjags)
library(rjags)
library(R2jags)
library(digest)
library(R2admb)
library(ggplot2)
library(pscl)
library(boot)
library(MASS)
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


#vetor de médias
mu<-exp(beta0+beta1*x1+beta*x2)
View(as.data.frame(mu))

#variável resposta
y<-NULL
for (i in 1:n) {
  y[i]<-rpois(1,mu[i])  
}
y

#descritiva
hist(y)
mean(y)
var(y)

#construindo o conj. de dados
data<-cbind(y,x1,x2)
data<-as.data.frame(data)

data
toc()

################################################ ANÁLISE FREQUENTISTA (cLÁSSICA)

tic("Ajuste utilizando o pacote base")

#base
p.glm<- glm(y ~ x1 +x2, family = poisson, data = data)
summary(p.glm)

toc()

################################################################################
#glmmTMB
#classico
if(!require(glmmTMB)) install.packages("glmmTMB");library(glmmTMB)

tic("Ajuste utilizando pacote glmmTMB") 
p.glmmTMB<-glmmTMB(y ~ x1 +x2, family=poisson, data=data)
summary(p.glmmTMB)
toc() 

################################################################################
#clasico
#mgcv
if(!require(mgcv)) install.packages("mgcv")
library(mgcv)
tic("Ajuste utilizando pacote mgcv")
p.mgcv = gam(y~x1+x2,family=poisson(),data=data)
summary(p.mgcv)
AIC(p.mgcv)
toc()

################################################################################

tic("Ajuste utilizando pacote glmmADMB") #bayesiano

#glmmADMB
p.glmmADMB<-glmmadmb(y ~ x1 +x2, family="poisson",data=data)
summary(p.glmmADMB)

toc()

################################################################################

tic("Ajuste utilizando pacote gamlss")

#gamlss
p.gamlss<-gamlss(y ~ x1 +x2, family=PO, data=data) 
summary(p.gamlss)

toc()


############################################################## ANÁLISE BAYESIANA

tic("Ajuste utilizando pacote arm")

#arm
p.arm<- bayesglm(y ~ x1 +x2, family = poisson, data = data)
summary(p.arm)
#BIC(p.arm)

toc()

################################################################################

tic("Ajuste utilizando o pacote jags")

#JAGS
datas <- list ("y"=y,"x1"=x1,"x2"=x2,"n"=n)  
inits <- function(){list(beta=c(0.5,1,-1))}
param<-c("beta","deviance")
p.jags<-run.jags(model="modelp.txt",monitor= param,data= datas,inits= inits,n.chains=3,
                 sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)

#print(p.jags)
summary(p.jags)
#summary(p.jags)[,4]
#extract(p.jags, what='dic')

toc()

################################################################################

tic("Ajuste utilizando pacote inla")

#inla
p.inla <- inla(y~x1+x2, family='poisson',
               data=data,control.family=list(link='log'),
               control.predictor=list(link=1, compute=TRUE),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(p.inla)
toc()

################################################################################

tic("Ajuste utilizando pacote MCMCglm")

#MCMCglmm
p.MCMCglmm<-MCMCglmm(y ~ x1 +x2, random=NULL, family="poisson",data=data, 
                     verbose=FALSE,nitt=4000, burnin=1000, thin=5)
summary(p.MCMCglmm)

toc()

################################################################################

tic("Ajuste utilizando pacote brms")

#brms
p.brm <- brm(y ~ x1 + x2, family = "poisson", data = data)
summary(p.brm)
WAIC(p.brm)

toc()

######################################################### Cálculo do afastamento

#CADA MODELO E DIFERENTE, DAR UMA OLHADA NA CHAMADA DOS COEFICIENTES

vetor_de_afastamento <- NULL

#FREQUENTISTA (cLÁSSICA)
(afastamento_glm <- (((0.5 - p.glm$coefficients[1])^2 + (1 - p.glm$coefficients[2])^2 + ((-1) - p.glm$coefficients[3])^2)/3))
vetor_de_afastamento[1] <- afastamento_glm
(afastamento_glmmADMB <- (((0.5 - p.glmmADMB$b[1])^2 + (1 - p.glmmADMB$b[2])^2 + ((-1) - p.glmmADMB$b[3])^2)/3))
vetor_de_afastamento[2] <- afastamento_glmmADMB
(afastamento_gamlss <- (((0.5 - p.gamlss$mu.coefficients[1])^2 + (1 - p.gamlss$mu.coefficients[2])^2 + ((-1) - p.gamlss$mu.coefficients[3])^2)/3))
vetor_de_afastamento[3] <- afastamento_gamlss
  
#ANÁLISE BAYESIANA
(afastamento_arm <- (((0.5 - p.arm$coefficients[1])^2 + (1 - p.arm$coefficients[2])^2 + ((-1) - p.arm$coefficients[3])^2)/3))
vetor_de_afastamento[4] <- afastamento_arm
(afastamento_jags <- (((0.5 - 0.4555522)^2 + (1 - 0.8338336)^2 + ((-1) - (-0.9656084))^2)/3)) 
vetor_de_afastamento[5] <- afastamento_jags
(afastamento_inla <- (((0.5 - 0.45800)^2 + (1 - 0.83200)^2 + ((-1) - (-0.96500))^2)/3))  
vetor_de_afastamento[6] <- afastamento_inla
(afastamento_MCMCglm <- (((0.5 - 0.50510)^2 + (1 - 0.70140)^2 + ((-1) - (-0.99100))^2)/3)) 
vetor_de_afastamento[7] <- afastamento_MCMCglm
(afastamento_brms <- (((0.5 - 0.46)^2 + (1 - 0.83)^2 + ((-1) - (-0.97))^2)/3))
vetor_de_afastamento[8] <- afastamento_MCMCglm
vetor_de_afastamento
