#### #### #### #### #### #### #### ####
#### Stima modelli Lineare - MCMC
#### #### #### #### #### #### #### ####
   

## ## ## ## ## ##
## Modello Lineare - MCMC
## ## ## ## ## ##


ModLinC = function(formula,  beta.mean, beta.variance,sigma2.a, sigma2.b, start.beta, start.sigma2, iter, burnin, thin, Data, sd.prop = 1)
{

  ### Definiamo l'osservazioni e la matrice X
  Cov     = model.matrix(formula,Data)
  NameObs = as.character(terms(formula)[[2]])
  Obs     = Data[,NameObs,drop=F] # drop serve a fare in modo che il risultato sia ancora una matrice

  nbeta = ncol(Cov)
  nobs  = nrow(Cov)

  ### definisco quanti campioni a posteriori devo salvare
  ### e gi oggetti che ritorna la funzione
  NsampleSave     = floor((iter-burnin)/thin)
  BetaSave        = matrix(NA,ncol=nbeta, nrow=NsampleSave )
  sigma2Save      = matrix(NA,ncol=1, nrow=NsampleSave )

  ### definisco i valori correnti di beta e sigma2
  betaMCMC   = start.beta
  sigma2MCMC = start.sigma2



  ## oggetti utili per l'MCMC
  appSamp = burnin
  i = 1
  XtX = t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
  Xty = t(Cov[i,,drop=F])*Obs[i,1]
  for(i in 2:nobs)
  {
    XtX = XtX+t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
    Xty = Xty+t(Cov[i,,drop=F])*Obs[i,1]
  }
  
  
  
  
  
  for(iMCMC in 1:NsampleSave)
  {
    for(jMCMC in 1:appSamp)
    {
      ## ## ## ## ## ## ## ##
      ## campioniamo Beta
      ## ## ## ## ## ## ## ##

      # prior parameters
      Vpost = solve(beta.variance)
      Mpost = solve(beta.variance)%*%matrix(beta.mean,ncol=1)

      # contributi delle osservazioni
      Vpost = Vpost+XtX/sigma2MCMC
      Mpost = Mpost+Xty/sigma2MCMC

      Vpost = solve(Vpost)
      Mpost = Vpost%*%Mpost

      #simlazione Gibbs
      betaMCMC = rmnorm(1,Mpost,Vpost)

      ## ## ## ## ## ## ## ##
      ## campioniamo sigma2 _ Metropolis
      ## ## ## ## ## ## ## ##

      # propongo  un nuovo valore
      sigma2prop = exp(rnorm(1,log(sigma2MCMC),sd.prop))

      # calcolo l'acceptance rate
      # potete scrivelo in due modi
      # 1) la distribuzione proposal è un log normale (esponenziale di una normale) e la prior gamma
      # 2) invece di sigma2 lavoriamo con tau=log(sigma2), che ha proposal simmetrica, e la prior la dobbiamo trovare
      # con la regola di trasformazioni di variabili:
      # f(tau) = f(sigma2(tau)) |d sigma2(tau) / d tau|
      # il second metodo, che è quello che usiamo, ha il vantaggio
      # che essendo la proposa simmetrica, si semplifica nel
      # calcolo di alpha

      # prior contribtion di tau
      logalpha = 0
      logalpha = logalpha+dgamma(sigma2prop, shape=sigma2.a, rate= sigma2.b,log=T)+log(sigma2prop)
      logalpha = logalpha-(dgamma(sigma2MCMC, shape=sigma2.a, rate= sigma2.b,log=T)+log(sigma2MCMC))

      # likelihood contribution
      # for(i in 1:nobs)
      # {
      #   logalpha = logalpha+dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2prop^0.5, log=T)
      #   logalpha = logalpha-dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2MCMC^0.5, log=T)
      # }

      logalpha = logalpha+sum(dnorm(Obs[1:nobs,],Cov[1:nobs,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2prop^0.5, log=T))
      logalpha = logalpha-sum(dnorm(Obs[1:nobs,],Cov[1:nobs,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2MCMC^0.5, log=T))
        
      alpha = min(1,exp(logalpha))

      u = runif(1,0,1)
      if(u<alpha)
      {
        sigma2MCMC = sigma2prop
      }
    }
    appSamp = thin

    BetaSave[iMCMC,]    = betaMCMC
    sigma2Save[iMCMC,]  = sigma2MCMC
  }
  return(list(Beta =BetaSave, sigma2 = sigma2Save ))
}


#### #### #### #### #### #### ####
#### TESTIAMO IL CODICE CON DATI SIMULATI
#### #### #### #### #### #### ####
set.seed(10)
n       = 100000
beta    = c(1,2,1,-1)
x1      = runif(n, 0,1)
x2      = runif(n, -1,1)
x3      = rnorm(n, 0,1)
x4      = rnorm(n, 0,1)
x5      = rnorm(n, 0,1)
x6      = rnorm(n, 0,1)
x7      = rnorm(n, 0,1)
x8      = rnorm(n, 0,1)
x9      = rnorm(n, 0,1)
x10      = rnorm(n, 0,1)
sigma2  = 1
y       = beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+rnorm(n, 0,sigma2^0.5)

Data = data.frame(y =y, x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9)

start_time <- Sys.time()
Results = ModLin(
  y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9,
  beta.mean = rep(0,10),
  beta.variance=diag(100,10),
  sigma2.a=1,
  sigma2.b=1,
  start.beta=rep(0,10),
  start.sigma2=1,
  iter=10000,
  burnin=1,
  thin=1,
  Data=Data
)
end_time <- Sys.time()

end_time-start_time

