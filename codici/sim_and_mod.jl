
#= Cose da fare la prima volta
La prima cosa da fare è creare il manifest e project 

using Pkg
cd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/julia/codici")
Pkg.instantiate() 

Pkg.instantiate va usato solo la prima volta che si crea il progetto
=#


#= Per il pacchetto

# funzione per creare il pacchetto

# PkgSkeleton.generate(pathtofile)

=#

#### #### #### #### #### 
#### Packages
#### #### #### #### #### 

using Pkg
Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/julia/codici")
# usare "]" e poi add Nome pacchetto, oppure
# Pkg.add(nome pacchetto)
# e
#  Pkg.develop(path=nome pacchetto) per il pacchetto che state sviluppando

# Pkg.add("Revise")
# Pkg.add("Distributions")
# Pkg.add("Random")
# Pkg.add("LinearAlgebra")
# Pkg.add("PDMats")
# Pkg.add("StatsBase")
# Pkg.develop(path="/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/julia/EsempioJulia/")

using Revise
using Distributions, Random
using LinearAlgebra, PDMats, StatsBase
using EsempioJulia
using RCall
using ToggleableAsserts
#### #### #### #### #### #### 
#### Data simulation
#### #### #### #### #### #### 

rngseed = 123456;
Random.seed!(rngseed);


n = 100::Int64
d = 10::Int64
ntot = Int64(n*d)
ncov_num  = 5::Int64
ncov_fact = 2::Int64
nfact = 3

# Parameters
sigma2 = 1.0::Float64

## creaiamo la matrice X

Xmat_app = Matrix{Float64}(undef, ntot,ncov_num+ ncov_fact )

for i = 1:ncov_num
    Xmat_app[:,i] = rand(Normal(0.0,2.0), ntot)
end

for i = (ncov_num+1):(ncov_num+ncov_fact)
    Xmat_app[:,i] = sample(1:nfact, ntot)
end

@rput ncov_num;
@rput ncov_fact;
@rput Xmat_app;

R"""

colnames(Xmat_app) = paste("x", 1:(ncov_num+ncov_fact), sep="")
Xmat_app = data.frame(Xmat_app)

for(i in (ncov_num+1):(ncov_num+ncov_fact))
{
    Xmat_app[,i] = factor(Xmat_app[,i])
}


form = "~"
for(i in 1:(ncov_num+ncov_fact))
{
    form = paste(form, " + x", i, sep="")
}

Xmat = model.matrix(formula(form), Xmat_app)


"""

@rget Xmat;

# se create ogetti "grandi", vi consiglio di elininarli da R, visto che la memoria non è condivisa
R"""

rm(Xmat)
rm(Xmat_app)

"""

## controllate sempre che le cose siano nella forma che volete
Xmat::Matrix{Float64};

## creaiamo uan funzione per generare i dati
# ma prima simuliamo i beta

regcoef = [rand(Normal(0.0,1.0)) for i = 1:size(Xmat,2)]::Vector{Float64};


ysim = sim(Xmat, regcoef,sigma2);

## oppure

ysim = Vector{Float64}(undef, ntot);
sim(ysim, Xmat, regcoef,sigma2);

# trasormiamo ysim in una matrice
ysim = reshape(ysim, (n,d))
### facciamo dei plot

@rput ysim;
@rput d
R"""

plot(ysim, col=rep(1:d, times=length(ysim)/d), pch=20 )

"""


#### #### #### #### #### 
#### Mutable and Immutable objects- examples
#### #### #### #### #### 

Mat_a = ones(Float64,5 ,5)
Mat_b = ones(Float64,5 ,5)

# check if they have the same values
Mat_a == Mat_b
# check if they are the same object
Mat_a === Mat_b

#
Mat_c = Mat_a

Mat_a == Mat_c
Mat_a === Mat_c

Mat_a[1,1] = 10.0
Mat_c[1,2] = -10.0
Mat_b[1,1:2]
Mat_a[1,1:2]
Mat_c[1,1:2]

Mat_d = deepcopy(Mat_a)
Mat_a[1,1] = 100.0

Mat_d[1,1]


#### #### #### #### #### 
#### view - examples
#### #### #### #### #### 

Mat = reshape(rand(Uniform(0.0,1.0),50^2),(50,50));

vect = Mat[1,[1,2,5]]
typeof(vect)
vect[3] = 10.0
Mat[1,5]

vect_view = @view  Mat[1,[1,2,5]]
typeof(vect_view)

vect_view[3] = 10.0
Mat[1,5]



## stimiamo il modello
toggle(false)

regcoefOUT, sigma2OUT = MCMC(
    ysim;
    X = Xmat,
    mcmc = (iter = 10000, thin = 2, burnin = 2),
    priors  = (regcoef=[0.0,1000.0], sigma2 = [0.1,0.1]),
    init  = (regcoef = [0.0 for i = 1:size(Xmat,2)], sigma2 = [1.0]),
    #sigma2Update = EsempioJulia.Sigma2Gibbs()
    sigma2Update = EsempioJulia.Sigma2Metropolis()
    # sigma2Update = EsempioJulia.Sigma2Hamil()
    #sigma2Update = EsempioJulia.Sigma2NoUp()

);

## test Packages
# Pkg.test("EsempioJulia")

@rput regcoefOUT;
@rput sigma2OUT;
@rput regcoef;
@rput sigma2;
R"""

pdf("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/julia/codici/MCMCout.pdf")

par(mfrow=c(3,3))
for(i in 1:nrow(regcoefOUT))
{
    
    plot(regcoefOUT[i,], type="l")
    abline(h =regcoefOUT[i], col=2, lwd=2)
    acf(regcoefOUT[i,])
    hist(regcoefOUT[i,])

}

plot(sigma2OUT[1,], type="l")
abline(h =sigma2, col=2, lwd=2)
acf(sigma2OUT[1,])
hist(sigma2OUT[1,])
dev.off()


"""










