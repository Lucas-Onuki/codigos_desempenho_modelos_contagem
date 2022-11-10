#Pacotes necessários:

if(!require(gamlss)) install.packages("gamlss")  #abordagem clássica
library(gamlss)

# remove.packages("StanHeaders")
# remove.packages("rstan")
# remove.packages("gamlss")
# remove.packages("brms")
# remove.packages("tictoc")
# file.remove(".RData")

# install.packages(c("StanHeaders","rstan"),type = "source")

if(!require(brms)) install.packages("brms") #abordagem bayesiana
library(brms)

# if(!require(rstan)) install.packages("rstan") #manipulação do tamanho (brms) 
library(rstan)
# library(StanHeaders)

if(!require(tictoc)) install.packages("tictoc") #para cálculo do tempo
library(tictoc)

# install.packages("doParallel") #pacote para rodar em pararelo
# library(doParallel)
# 
# install.packages("foreach") 
# library(foreach)

#Preparando para rodar em pararelo
# numCores <- detectCores()
# numCores
# 
# registerDoParallel(numCores)

###############################################################################
# Função para simular um conjunto de dados a partir do modelo ZINB e          #
# Ajuste usando modelo ZINB a paritir de gamlss (Clássico) e BRMS (bayesiano) #
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
  
  
  #gerando dados para glm (---PARÂMETROS DE EXEMPLO---)
  # n <- 500
  # beta <- c(0.5,1,-1)
  # delta <- c(-4,0.5,-3.5)
  
  data<-geraZINB(n,beta,delta) 
  dataLglm=data.frame(data)
  names(dataLglm)<-c("y","x1","x2")
  #dataLglm
  
  ### AJUSTES ----
  
  ##GAMLSS
  
  tic("Ajuste utilizando pacote gamlss")
  zinb.gamlss<-gamlss(formula=y ~ x1 +x2, nu.formula =y ~ x1 +x2, family = ZINBI, data = dataLglm,
                      trace =TRUE)
  summary(zinb.gamlss)
  t1 = toc()
  
  #tempo
  tempo_glm = t1[[2]] - t1[[1]]
  
  #size estimado
  sigma<-exp(zinb.gamlss$sigma.coefficients)
  sigma_output <- 1/sigma
  
  #coeficientes de mu + size
  gamlssfitglm <- c(zinb.gamlss$mu.coefficients,zinb.gamlss$nu.coefficients)
  gamlssfitglm
  
  ##BRMS
  
  tic("Ajuste utilizando pacote brms")
  zinb.brm <- brm(bf(y ~ x1 + x2, zi ~ x1+x2), 
                  data = dataLglm, family = zero_inflated_negbinomial(),
                  prior=set_prior("normal(0,1)"), seed=170400963, refresh=500)
  summary(zinb.brm)
  t2 = toc()
  
  #coeficientes de mu
  brmfitglm <- fixef(zinb.brm)
  n_brmfitglm <- brmfitglm[,1]
  
  #tempo
  tempo_brm = t2[[2]] - t2[[1]]
  
  #size
  list <- extract(zinb.brm$fit)
  tamanho_estimado_brm <- mean(list$shape)
  
  #Output final
  output<-c(n,sigma_output,tempo_glm,gamlssfitglm,tamanho_estimado_brm,tempo_brm,n_brmfitglm)
  
  return(output)
}

#EXEMPLO----
# n <- 500
# beta <- c(0.5,1,-1)
# delta <- c(-4,0.5,-3.5)
# 
# modeloZINB_two_parts(n,beta,delta)

# 3betas
# 3deltas
# tamanho e tempo dos métodos

###################################################################################
# Função para simular conjunto de dados com réplicas a partir de um modelo ZINB e #
# Ajuste usando modelo ZINB a paritir de gamlss (Clássico) e BRMS (bayesiano)     #
###################################################################################


replicZINB_two_parts<-function(M,n,beta,delta) #M sendo o número de réplicas
{
  replica = matrix(0,M,17) #matriz que armazena os resultados
  
  colnames(replica) <- c("size",
                       "size_estimado_glm",
                       "tempo_glm",
                       "Beta_0glm",
                       "Beta_1glm",
                       "Beta_2glm",
                       "Delta_0glm",
                       "Delta_1glm",
                       "Delta_2glm",
                       "size_estimado_glm",
                       "tempo_brm",
                       "Beta_0brm",
                       "Delta_0brm",
                       "Beta_1brm",
                       "Beta_2brm",
                       "Delta_1brm",
                       "Delta_2brm")
  
  for(i in 1:M){
    replica[i,] = modeloZINB_two_parts(n,beta,delta) #ajustes
  }
  
  #médias ----
  
  #Gamlss ----
  
  mean_B0_glm = mean(replica[,4])
  mean_B1_glm = mean(replica[,5])
  mean_B2_glm = mean(replica[,6])
  
  mean_D0_glm = mean(replica[,7])
  mean_D1_glm = mean(replica[,8])
  mean_D2_glm = mean(replica[,9])
  
  #tempo computacional
  media_do_tempo_glm = mean(replica[,3])
  
  #size
  media_tamanho_glm = mean(replica[,2]) 
  
  ###############################-
  
  #Brms ----
  
  mean_B0_brm = mean(replica[,12])
  mean_B1_brm = mean(replica[,14])
  mean_B2_brm = mean(replica[,15])
  
  mean_D0_brm = mean(replica[,13])
  mean_D1_brm = mean(replica[,16])
  mean_D2_brm = mean(replica[,17])
  
  #tempo computacional
  media_do_tempo_brm = mean(replica[,11])
  
  #size
  media_tamanho_brm = mean(replica[,10])
  
  ###############################-
  
  #Viés ----
  
  #Gamlss ----
  
  bias_B0_glm = mean_B0_glm - beta[1]
  bias_B1_glm = mean_B1_glm - beta[2]
  bias_B2_glm = mean_B2_glm - beta[3]
   
  bias_D0_glm = mean_D0_glm - delta[1]
  bias_D1_glm = mean_D1_glm - delta[2]
  bias_D2_glm = mean_D2_glm - delta[3]
  
  ###############################-
  
  #Brms ----
  
  bias_B0_brm = mean_B0_brm - beta[1]
  bias_B1_brm = mean_B1_brm - beta[2]
  bias_B2_brm = mean_B2_brm - beta[3]
   
  bias_D0_brm = mean_D0_brm - delta[1]
  bias_D1_brm = mean_D1_brm - delta[2]
  bias_D2_brm = mean_D2_brm - delta[3]
  
  ###############################-
  
  #Desvio Padrão ----
  
  #Gamlss ----
  
  sd_B0_glm = sd(replica[,4])
  sd_B1_glm = sd(replica[,5])
  sd_B2_glm = sd(replica[,6])
  
  sd_D0_glm = sd(replica[,7])
  sd_D1_glm = sd(replica[,8])
  sd_D2_glm = sd(replica[,9])
  
  #tempo computacional
  sd_tempo_glm = sd(replica[,3])
  
  #size
  sd_tamanho_glm = sd(replica[,2])
  
  ###############################-
  
  #Brms ----
  
  sd_B0_brm = sd(replica[,12])
  sd_B1_brm = sd(replica[,14])
  sd_B2_brm = sd(replica[,15])
  
  sd_D0_brm = sd(replica[,13])
  sd_D1_brm = sd(replica[,16])
  sd_D2_brm = sd(replica[,17])

  #tempo computacional
  sd_tempo_brm = sd(replica[,11])
  
  #size
  sd_tamanho_brm = sd(replica[,10])
  
  ###############################-
  
  #MSE ----
  
  #Gamlss ----
  
  mse_B0_glm = mean((replica[,2]-beta[1])^2)
  mse_B1_glm = mean((replica[,3]-beta[2])^2)
  mse_B2_glm = mean((replica[,4]-beta[3])^2)
  
  mse_D0_glm = mean((replica[,5]-delta[1])^2)
  mse_D1_glm = mean((replica[,6]-delta[2])^2)
  mse_D2_glm = mean((replica[,7]-delta[3])^2)
  
  ###############################-
  
  #Brms ----
  
  mse_B0_brm = mean((replica[,8]-beta[1])^2)
  mse_B1_brm = mean((replica[,9]-beta[2])^2)
  mse_B2_brm = mean((replica[,10]-beta[2])^2)
  
  mse_D0_brm = mean((replica[,11]-delta[2])^2)
  mse_D1_brm = mean((replica[,12]-delta[2])^2)
  mse_D2_brm = mean((replica[,13]-delta[2])^2)
  
  ###############################-
  
  resmean = c(media_tamanho_glm,
              media_tamanho_brm,
              media_do_tempo_glm,
              media_do_tempo_brm,
              mean_B0_glm,
              mean_B1_glm,
              mean_B2_glm,
              mean_D0_glm,
              mean_D1_glm,
              mean_D2_glm,
              mean_B0_brm,
              mean_B1_brm,
              mean_B2_brm,
              mean_D0_brm,
              mean_D1_brm,
              mean_D2_brm)
  
  names(resmean)<-c("media_tamanho_glm",
                    "media_tamanho_brm",
                    "media_do_tempo_glm",
                    "media_do_tempo_brm",
                    "mean_B0_glm",
                    "mean_B1_glm",
                    "mean_B2_glm",
                    "mean_D0_glm",
                    "mean_D1_glm",
                    "mean_D2_glm",
                    "mean_B0_brm",
                    "mean_B1_brm",
                    "mean_B2_brm",
                    "mean_D0_brm",
                    "mean_D1_brm",
                    "mean_D2_brm")
  
  resbias = c(bias_B0_glm,
              bias_B1_glm,
              bias_B2_glm,
              bias_D0_glm,
              bias_D1_glm,
              bias_D2_glm,
              bias_B0_brm,
              bias_B1_brm,
              bias_B2_brm,
              bias_D0_brm,
              bias_D1_brm,
              bias_D2_brm)
  
  names(resbias)<-c("bias_B0_glm",
                    "bias_B1_glm",
                    "bias_B2_glm",
                    "bias_D0_glm",
                    "bias_D1_glm",
                    "bias_D2_glm",
                    "bias_B0_brm",
                    "bias_B1_brm",
                    "bias_B2_brm",
                    "bias_D0_brm",
                    "bias_D1_brm",
                    "bias_D2_brm")
  
  ressd = c(sd_tamanho_glm,
            sd_tamanho_brm,
            sd_tempo_glm,
            sd_tempo_brm,
            sd_B0_glm,
            sd_B1_glm,
            sd_B2_glm,
            sd_D0_glm,
            sd_D1_glm,
            sd_D2_glm,
            sd_B0_brm,
            sd_B1_brm,
            sd_B2_brm,
            sd_D0_brm,
            sd_D1_brm,
            sd_D2_brm)
  
  names(ressd)<-c("sd_tamanho_glm",
                  "sd_tamanho_brm",
                  "sd_tempo_glm",
                  "sd_tempo_brm",
                  "sd_B0_glm",
                  "sd_B1_glm",
                  "sd_B2_glm",
                  "sd_D0_glm",
                  "sd_D1_glm",
                  "sd_D2_glm",
                  "sd_B0_brm",
                  "sd_B1_brm",
                  "sd_B2_brm",
                  "sd_D0_brm",
                  "sd_D1_brm",
                  "sd_D2_brm")
  
  resmse = c(mse_B0_glm,
             mse_B1_glm,
             mse_B2_glm,
             mse_D0_glm,
             mse_D1_glm,
             mse_D2_glm,
             mse_B0_brm,
             mse_B1_brm,
             mse_B2_brm,
             mse_D0_brm,
             mse_D1_brm,
             mse_D2_brm)
  
  names(resmse) <- c("mse_B0_glm",
                     "mse_B1_glm",
                     "mse_B2_glm",
                     "mse_D0_glm",
                     "mse_D1_glm",
                     "mse_D2_glm",
                     "mse_B0_brm",
                     "mse_B1_brm",
                     "mse_B2_brm",
                     "mse_D0_brm",
                     "mse_D1_brm",
                     "mse_D2_brm")
  
  b <- list(repl = replica, 
            mean = resmean, 
            bias = resbias, 
            sd = ressd, 
            mse = resmse)
  print(b)
}


M = 50
n = 1000


# tic("replicas")
# replicZINB_two_parts(M,n,beta,delta)
# ttt = toc()

#Para rodar de uma vez ----

a <- list()

#tamanhos amostrais
n <- c(200,500,1000)

#parâmetros -
beta <- c(0.5,1,-1) 
delta <- c(-4,0.5,-3.5)

for (i in 1:length(n)) {
  a[[i]] <- replicZINB_two_parts(M,n[i],beta,delta)
}

#para encerrar
# stopImplicitCluster()
