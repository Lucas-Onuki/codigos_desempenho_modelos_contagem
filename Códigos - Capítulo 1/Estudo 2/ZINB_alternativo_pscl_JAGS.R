#Pacotes necessários:

if(!require(tictoc)) install.packages("tictoc") #para cálculo do tempo
library(tictoc)

###############################################################################
# Função para simular um conjunto de dados a partir do modelo ZINB e          #
# Ajuste usando modelo ZINB a paritir de pscl (Clássico) e JAGS (bayesiano)   #
###############################################################################

modeloZINB_two_parts<- function(n,beta,delta){
  
  #função para gerar dados zero-inflacionados segundo binomial negativa
  geraZINB<-function(n,beta,delta)
  {
    n <- n #tamanho amostral
    
    #parâmetros 
    #betas (parte binomial negativa)----
    beta2 <- beta[3]
    beta1 <- beta[2] 
    beta0 <- beta[1]
    #deltas (parte inflacionada)----
    delta2 <- delta[3]
    delta1 <- delta[2] 
    delta0 <- delta[1]
    
    x <- matrix(rnorm(2*n,0,1),ncol=2) #geração das covariáveis
    mu <- exp(beta0 + x%*%beta[-1])    #vetor de médias
    pL <- plogis(delta0+x%*%delta[-1]) #probabilidade de zeros
    
    #geração da variável resposta----
    y <- NULL
    for (i in 1:n) {
      y.zero <- rbinom(n, prob =  1-pL[i], size = 1)
      y[i] <- ifelse(y.zero == 0, 0, rnbinom(n = 1, mu = mu[i], size = 2)) #phi = size (tamanho)
    }
    
    return(cbind(y,x))
  }
  
  
  #gerando dados para pscl (---PARÂMETROS DE EXEMPLO---)
  n <- 500
  beta <- c(0.5,1,-1)
  delta <- c(-4,0.5,-3.5)
  
  data<-geraZINB(n,beta,delta) 
  dataLpscl=data.frame(data)
  names(dataLpscl)<-c("y","x1","x2")
  #dataLpscl
  
  ### AJUSTES ----
  ###############################################################################
  
  ##PSCL
  if(!require(pscl)) install.packages("pscl")  #abordagem clássica
  library(pscl) 
  
  tic("Ajuste utilizando o pacote pscl")
  zinb.zeroinfl <- zeroinfl(y ~ x1 + x2 | x1 + x2,  data = dataLpscl, dist = "negbin")
  summary(zinb.zeroinfl)
  t1 = toc()
  
  zinb.zeroinfl$coefficients
  #size estimado
  zinb.zeroinfl$theta
  
  #tempo
  tempo_pscl = t1[[2]] - t1[[1]]
  
  #parametros
  psclfitpar <- c(zinb.zeroinfl$coefficients$count,zinb.zeroinfl$coefficients$zero,zinb.zeroinfl$theta)
  psclfitpar
  
  ###############################################################################
  
  ##JAGS
  
  if(!require(runjags)) install.packages("runjags") #abordagem bayesiana
  library(runjags)
  
 # library("R2jags")
 # setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
# Model in the JAGS format
  modelzinb<- "model {
  	for(i in 1:n){
	       y[i] ~ dnegbin(p[i],phi)
               p[i] <- phi/(phi+(1-zero[i])*lambda[i])- 1e-10*zero[i]
              	#lambda[i]<- max(0.00000000001,lambdaaux[i])
		lambda[i]<- exp(beta[1]+beta[2]*x1[i]+beta[3]*x2[i])

              ## Zero-Inflation
   		zero[i]~dbern(q[i])
		q[i]<- max(0.00000000001,min(0.99999999999,1-qaux[i]))
		logit(qaux[i])<- alpha[1]+alpha[2]*x1[i] +alpha[3]*x2[i]
		
            
	}

#priors
phi~dgamma(0.1,0.01)
for (j in 1:3){
	beta[j]~dnorm(0,0.01)
	alpha[j]~dnorm(0,0.01)
	}
}"
  

  attach(dataLpscl)
  datas <- list ("y"= y,"x1"=x1,"x2"=x2,"n"=n)  
  inits <- function(){list(beta=c(0.5,1,-1),alpha=c(-4,0.5,-3.5), phi=2)} 
  param<-c("beta","alpha","phi")
  
  tic("Ajuste utilizando o pacote jags")
  zinb.jags2<- run.jags(model=modelzinb ,monitor= param,data= datas,inits= inits,n.chains=3,
                             sample=1000,adapt=3000,burnin=1000,thin=5,modules='runjags',jags.refresh=40)
  summary(zinb.jags2)
  t2 =  toc()
  #tempo
  tempo_jags = t2[[2]] - t2[[1]] 
  
  
  library('coda')
  mcmclist <- as.mcmc.list(zinb.jags2, vars=c("beta","alpha","phi"))
  
  if(!require(MCMCvis)) install.packages("MCMCvis")
  library(MCMCvis)
  resfitpar<-MCMCsummary(mcmclist, round = 7)
  resfitpar
  jagsfitpar<-resfitpar$mean
  jagsfitpar
  
  #Output final
  output<-c(n,tempo_pscl,psclfitpar,tempo_jags,jagsfitpar)
  
  return(output)
}

replicZINB_two_parts<-function(M,n,beta,delta)
{
  replica=matrix(0,M,17)
  colnames(replica)<-c("size","tempo_pscl","Beta_0pscl","Beta_1pscl","Beta_2pscl","Delta_0pscl","Delta_1pscl","Delta_2pscl","size_estimado_pscl","tempo_jags","Beta_0jags","Delta_0jags","Beta_1jags","Beta_2jags","Delta_1jags","Delta_2jags","size_estimado_pscl")
  for(i in 1:M){
    replica[i,]=modeloZINB_two_parts(n,beta,delta)
  }
  
  #médias ----
  mean_B0_pscl=mean(replica[,3])
  mean_B1_pscl=mean(replica[,4])
  mean_B2_pscl=mean(replica[,5])
  
  mean_D0_pscl=mean(replica[,6])
  mean_D1_pscl=mean(replica[,7])
  mean_D2_pscl=mean(replica[,8])
  
  mean_B0_jags=mean(replica[,11])
  mean_B1_jags=mean(replica[,12])
  mean_B2_jags=mean(replica[,13])
  
  mean_D0_jags=mean(replica[,14])
  mean_D1_jags=mean(replica[,15])
  mean_D2_jags=mean(replica[,16])
  
  media_do_tempo_pscl=mean(replica[,2])
  media_do_tempo_jags=mean(replica[,10])
  
  media_tamanho_pscl=mean(replica[,9])  
  media_tamanho_jags=mean(replica[,17])
  
  
  #Viés ----
  bias_B0_pscl=mean_B0_pscl-beta[1]
  bias_B1_pscl=mean_B1_pscl-beta[2]
  bias_B2_pscl=mean_B2_pscl-beta[3]
  
  bias_B0_jags=mean_B0_jags-beta[1]
  bias_B1_jags=mean_B1_jags-beta[2]
  bias_B2_jags=mean_B2_jags-beta[3]
  
  bias_D0_pscl=mean_D0_pscl-delta[1]
  bias_D1_pscl=mean_D1_pscl-delta[2]
  bias_D2_pscl=mean_D2_pscl-delta[3]
  
  bias_D0_jags=mean_D0_jags-delta[1]
  bias_D1_jags=mean_D1_jags-delta[2]
  bias_D2_jags=mean_D2_jags-delta[3]
  
  #Desvio Padrão ----
  
  sd_tamanho_pscl=sd(replica[,9])/sqrt(M)
  sd_tamanho_jags=sd(replica[,17])/sqrt(M)
  
  sd_B0_pscl=sd(replica[,3])/sqrt(M)
  sd_B1_pscl=sd(replica[,4])/sqrt(M)
  sd_B2_pscl=sd(replica[,5])/sqrt(M)
  
  sd_D0_pscl=sd(replica[,6])/sqrt(M)
  sd_D1_pscl=sd(replica[,7])/sqrt(M)
  sd_D2_pscl=sd(replica[,8])/sqrt(M)
  
  sd_B0_jags=sd(replica[,11])/sqrt(M)
  sd_B1_jags=sd(replica[,12])/sqrt(M)
  sd_B2_jags=sd(replica[,13])/sqrt(M)
  
  sd_D0_jags=sd(replica[,14])/sqrt(M)
  sd_D1_jags=sd(replica[,15])/sqrt(M)
  sd_D2_jags=sd(replica[,16])/sqrt(M)
  
  sd_tempo_pscl=sd(replica[,2])/sqrt(M)
  sd_tempo_jags=sd(replica[,10])/sqrt(M)
  
  #MSE ----
  mse_B0_pscl<-mean((replica[,3]-beta[1])^2)
  mse_B1_pscl<-mean((replica[,4]-beta[2])^2)
  mse_B2_pscl<-mean((replica[,5]-beta[3])^2)
  
  mse_D0_pscl<-mean((replica[,6]-delta[1])^2)
  mse_D1_pscl<-mean((replica[,7]-delta[2])^2)
  mse_D2_pscl<-mean((replica[,8]-delta[3])^2)
  
  mse_B0_jags<-mean((replica[,11]-beta[1])^2)
  mse_B1_jags<-mean((replica[,12]-beta[2])^2)
  mse_B2_jags<-mean((replica[,13]-beta[2])^2)
  
  mse_D0_jags<-mean((replica[,14]-delta[2])^2)
  mse_D1_jags<-mean((replica[,15]-delta[2])^2)
  mse_D2_jags<-mean((replica[,16]-delta[2])^2)
  
  
  resmean=c(media_tamanho_pscl,media_tamanho_jags,media_do_tempo_pscl,media_do_tempo_jags,mean_B0_pscl,mean_B1_pscl,mean_B2_pscl,mean_D0_pscl,mean_D1_pscl,mean_D2_pscl,mean_B0_jags,mean_B1_jags,mean_B2_jags,mean_D0_jags,mean_D1_jags,mean_D2_jags)
  names(resmean)<-c("media_tamanho_pscl","media_tamanho_jags","media_do_tempo_pscl","media_tempo_jags","mean_B1_pscl","mean_B2_pscl","mean_D0_pscl","mean_D1_pscl","mean_D2_pscl","media_do_tempo_jags","mean_B0_jags","mean_B1_jags","mean_B2_jags","mean_D0_jags","mean_D1_jags","mean_D2_jags")
  
  resbias=c(bias_B0_pscl,bias_B1_pscl,bias_B2_pscl,bias_D0_pscl,bias_D1_pscl,bias_D2_pscl,bias_B0_jags,bias_B1_jags,bias_B2_jags,bias_D0_jags,bias_D1_jags,bias_D2_jags)
  names(resbias)<-c("bias_B0_pscl","bias_B1_pscl","bias_B2_pscl","bias_D0_pscl","bias_D1_pscl","bias_D2_pscl","bias_B0_jags","bias_B1_jags","bias_B2_jags","bias_D0_jags","bias_D1_jags","bias_D2_jags")
  
  ressd=c(sd_tamanho_pscl,sd_tamanho_jags,sd_tempo_pscl,sd_tempo_jags,sd_B0_pscl,sd_B1_pscl,sd_B2_pscl,sd_D0_pscl,sd_D1_pscl,sd_D2_pscl,sd_B0_jags,sd_B1_jags,sd_B2_jags,sd_D0_jags,sd_D1_jags,sd_D2_jags)
  names(ressd)<-c("sd_tamanho_pscl","sd_tamanho_jags","sd_tempo_pscl","sd_tempo_jags","sd_B0_pscl","sd_B1_pscl","sd_B2_pscl","sd_D0_pscl","sd_D1_pscl","sd_D2_pscl","sd_B0_jags","sd_B1_jags","sd_B2_jags","sd_D0_jags","sd_D1_jags","sd_D2_jags")
  
  resmse=c(mse_B0_pscl,mse_B1_pscl,mse_B2_pscl,mse_D0_pscl,mse_D1_pscl,mse_D2_pscl,mse_B0_jags,mse_B1_jags,mse_B2_jags,mse_D0_jags,mse_D1_jags,mse_D2_jags)
  names(resmse)<-c("mse_B0_pscl","mse_B1_pscl","mse_B2_pscl","mse_D0_pscl","mse_D1_pscl","mse_D2_pscl","mse_B0_jags","mse_B1_jags","mse_B2_jags","mse_D0_jags","mse_D1_jags","mse_D2_jags")
  
  b<-list(repl=replica,mean=resmean,bias=resbias,sd=ressd,mse=resmse)
  print(b)
}

M=50
n=1000
beta <- c(0.5,1,-1)
delta <- c(-4,0.5,-3.5)

tic("replicas")
replicZINB_two_parts(M,n,beta,delta)
ttt = toc()

#500 - 1000 replicas
#200 - 1000 replicas

#Em teoria com menos amostra, bayesiano é melhor


#Para rodar de uma vez ----
# aux = c(1,2,3,4,5,6,7,8,9,10)
# a <- NULL
# b <- NULL
# 
# for (i in aux) {
#   a[i] <- replicZINB_two_parts(M,200,beta,delta)
#   b[i] <- replicZINB_two_parts(M,500,beta,delta)
# }
