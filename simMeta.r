#simMeta
#Written by R. Graham Reynolds and Ben M. Fitzpatrick 2010
#www.rgrahamreynolds.info
#Citation: Reynolds, R.G. (2011) Islands, metapopulations, and archipelagos: genetic equilibrium and non-equilibrium dynamics of structured populations in the context of conservation. Ph.D. Dissertation, University of Tennessee. 


Step 1) Function File
R file for simulation of population data 

This file sets up the functions to simulate population models from Chapter I
Functions: LG; SMM; pop; smd; migrate; Fst, MMD
Step 1) Execute this file
Step 2) Add a Population Model File: either MM (Bottleneck or Extinction and 	Recolonization) or FDM (Bottleneck or Extinction/Recolonization)
Step 3) Add and adjust the migration matrix
Step 4) Add Simulation File- input parametrs and # of simulations
Execute all of the above

Notations are in UTK Orange
R code is in black
Functions are in light blue

# LG is the logistic growth function
LG=function(N,r,K,T){
	plot(0,0,type="n",xlim=c(0,T),ylim=c(0,2*K),xlab="generation",ylab="N")
	for(i in 1:T){
		N=N+r*N*(K-N)
		points(i,N)
		print(N)
	}}

#Function SMM simulates the stepwise mutation model for SSR data
SMM=function(n,N,theta,T,P){
	mu=theta/(4*N)
	for(i in 1:T){
		mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
		P=sample(P,2*N,TRUE)
		P=P+mutations
		}
	cbind(sample(P,n),sample(P,n))
	}

# make genepop file
pop=function(L=20,n=20,N=1000,theta=4,T=10000,P=rep(50,2000)){
	X=SMM(n,N,theta,T,P)
	for(i in 2:L){
		X=cbind(X,SMM(n,N,theta,T,P))}
	X
	}

#Function smd simulates a stepping stone metapopulation
smd=function(N,mu,P){
	#one generation of stepwise mutation (no drift)
	mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
	P+mutations
	}

#Function migrate is for one generation of migration among pops according to migration matrix m
migrate=function(pops,N,m){
	#pops is a matrix of k populations and L allele frequencies
	#N causes each population's contribution of migrants to be proportional to its size
	#npops is then the expected number of gametes with each allele in each population
	k=nrow(pops)
	npops=pops
	for(i in 1:k){
		npops[i,]=t(m[i,]*N)%*%pops
		}
	npops
	}

# Function Fst calculates pairwise Fst with no bias correction
Fst=function(x){
	# assume pops is a matrix with population membership in the first column and two columns of alleles
	tab=table(pop=rep(x[,1],2),allele=c(x[,2],x[,3]))
	k=nrow(tab)
	n=rowSums(tab)
	tab=tab/n
	f=matrix(NA,k,k)
	for(i in 1:k){
		for(j in 1:k){
		ptab=tab[c(i,j),]
		Ht=1-sum(colMeans(ptab)^2)
		Hs=1-sum(colMeans(ptab^2))
		f[i,j]=(Ht-Hs)/Ht
		}
		}
	f
	}

#Function MMD simulates T generations of mutation, migration, and drift in a set of populations
#this version works for exactly 12 populations
#pops is a table of allele frequencies
#N is a vector of initial population sizes
#r is the intrinsic growth rate
#K is a vector of carrying capacities
MMD=function(pops,N,m,mu,T,r,K,n){
	k=nrow(pops)
	alleles=colnames(pops)
	#plot(rep(0,k),N,xlim=c(0,T),ylim=c(0,2*max(K)),type="n",xlab="generation",ylab="N")
	for(i in 1:T){
		#migration
		popsm=migrate(pops,N,m)
		N=rowSums(popsm)
		#drift+mutation
		N=N+r*N*(K-N)
		N=round(N)
		N=replace(N,N<1,1)
		#points(rep(i,k),N)
		P=list(rep(NA,2*N[1]),rep(NA,2*N[2]),rep(NA,2*N[3]),rep(NA,2*N[4]),rep(NA,2*N[5]),rep(NA,2*N[6]),rep(NA,2*N[7]),rep(NA,2*N[8]),rep(NA,2*N[9]),rep(NA,2*N[10]),rep(NA,2*N[11]),rep(NA,2*N[12]))
		pnames=list(rep(1,2*N[1]),rep(2,2*N[2]),rep(3,2*N[3]),rep(4,2*N[4]),rep(5,2*N[5]),rep(6,2*N[6]),rep(7,2*N[7]),rep(8,2*N[8]),rep(9,2*N[9]),rep(10,2*N[10]),rep(11,2*N[11]),rep(12,2*N[12]))
		for(j in 1:k){
			P[[j]]=smd(N[j],mu,as.numeric(sample(alleles,size=2*N[j],replace=TRUE,prob=popsm[j,])))
			}
		alleles=names(table(unlist(P)))
		pops=table(unlist(pnames),unlist(P))/(2*N)	
		}
	n=rep(n,k)
	for(g in 1:12){if(n[g]>N[g]) n[g]=N[g]}
	p=list(rep(NA,2*n[1]),rep(NA,2*n[2]),rep(NA,2*n[3]),rep(NA,2*n[4]),rep(NA,2*n[5]),rep(NA,2*n[6]),rep(NA,2*n[7]),rep(NA,2*n[8]),rep(NA,2*n[9]),rep(NA,2*n[10]),rep(NA,2*n[11]),rep(NA,2*n[12]))
	pnames=list(rep(1,n[1]),rep(2,n[2]),rep(3,n[3]),rep(4,n[4]),rep(5,n[5]),rep(6,n[6]),rep(7,n[7]),rep(8,n[8]),rep(9,n[9]),rep(10,n[10]),rep(11,n[11]),rep(12,n[12]))
	for(j in 1:k){
			p[[j]]=sample(alleles,size=2*n[j],replace=TRUE,prob=pops[j,])
			}
	sam=matrix(as.numeric(unlist(p)),ncol=2,byrow=T)
	f=Fst(cbind(unlist(pnames),sam))
	list(pops=pops, sam=sam, Fst=f, N=N,pnames=unlist(pnames))
		}


#Step 2) Add a Population Model File: either MM (Bottleneck, Extinction/Recolonization, Equilibrium), or FDM (Bottleneck, Extinction/Recolonization, Varying Ne).
#Step 3) Add and adjust the Migration Matrix File or use the no migration matrix
#Step 4) Add Simulation File- input parameters and # of simulations
#Execute all of the above






Step 2) Population Model Files
Focal Deme Model- Bottleneck
Focal Deme Model- Extinction - Recolonization
Focal Deme Model- Varying Ne
Metapopulation Model (MM)- Equilibrium
Metapopulation Model (MM)- Bottleneck
Metapopulation Model (MM)- Extinction/Recolonization

I. Focal Deme Model- Bottleneck
# function for focal deme model with bottleneck
FDMbot=function(r,K,L,m,mu,time,n,b){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		bpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 drops to b individuals
		N2=c(b,rep(K,11))
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r,K=K,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

II. Focal Deme Model- Extinction - Recolonization
# function for focal deme model with extinction - recoloniztion
FDMer=function(r,K,L,m,mu,time,n){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		eqpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 goes extinct
		bpop=eqpop
		bpop[1,]=0
		N2=c(1,rep(K,11))
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r,K=K,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}


III. Focal Deme Model with Varying Ne 
# function for focal deme model with different Ne for focal deme 
FPMN=function(r,K=Ne,N=Nb,L,m,mu,time,n,b=Nb){
	pops=table(1:12,rep(50,12))
	N1=c(N,rep(K,11))
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		bpop=MMDc(pops=pops,N=N1,m=m,mu=mu,T=10*K,n=n)$pops
		# step 2: pop1 (Focal Deme) drops to b individuals
		N2=c(b,rep(K,11))
		#step 3: simulate more generations
		npop=MMDc(pops=bpop,N=N2,m=m,mu=mu,T=time,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

IV. Metapopulation Model- Equilibrium
#function for Metapopulation Model with constant population size (no demographic event)
MMeq=function(pops,N,m,mu,T,n){
	#simulate T generations of mutation, migration, and drift in a set of populations
	#this version works for exactly 12 populations
	#pops is a table of allele frequencies
	#N is a vector of initial population sizes
	k=nrow(pops)
	alleles=colnames(pops)
	for(i in 1:T){
		#migration
		popsm=migrate(pops,N,m)
		#drift+mutation
		P=list(rep(NA,2*N[1]),rep(NA,2*N[2]),rep(NA,2*N[3]),rep(NA,2*N[4]),rep(NA,2*N[5]),rep(NA,2*N[6]),rep(NA,2*N[7]),rep(NA,2*N[8]),rep(NA,2*N[9]),rep(NA,2*N[10]),rep(NA,2*N[11]),rep(NA,2*N[12]))
		pnames=list(rep(1,2*N[1]),rep(2,2*N[2]),rep(3,2*N[3]),rep(4,2*N[4]),rep(5,2*N[5]),rep(6,2*N[6]),rep(7,2*N[7]),rep(8,2*N[8]),rep(9,2*N[9]),rep(10,2*N[10]),rep(11,2*N[11]),rep(12,2*N[12]))
		for(j in 1:k){
	P[[j]]=smd(N[j],mu,as.numeric(sample(alleles,size=2*N[j],replace=TRUE,prob=popsm[j,])))
			}
		alleles=names(table(unlist(P)))
		pops=table(unlist(pnames),unlist(P))/(2*N)	
		}
	
	p=list(rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n))
	pnames=list(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n),rep(11,n),rep(12,n))
	for(j in 1:k){
			p[[j]]=sample(alleles,size=2*n,replace=T,prob=pops[j,])
			}
	sam=matrix(as.numeric(unlist(p)),ncol=2,byrow=T)
	f=Fst(cbind(unlist(pnames),sam))
	list(pops=pops, sam=sam, Fst=f, N=N, pnames=unlist(pnames))
		}


V. Metapopulation Model- Bottleneck
# Function MPMb is for a Metapopulation Model (MM) with bottlenecks
MPMbot=function(r,K,L,m,mu,Time,n,b){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	f=matrix(0,ncol=12,nrow=12)
	events=length(Time)
	island=NA
	# first locus
	# step 1: simulate for 10N generations to achieve equilibrium
	npop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)
	# step 2: set population to bottleneck size
	for(j in 1:events){
		EX=j
		island[j]=EX
		bpop=npop$pops
		N2=npop$N
		N2[EX]=b
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=Time[j],r=r,K=K,n=n)
		}
	Exdat=cbind(npop$pnames,npop$sam)
	f=f+Fst(cbind(npop$pnames,npop$sam))
	for(i in 2:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		npop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)
		# step 2: a random popultion goes extinct at each Time
		for(j in 1:events){
			EX=j
			island[j]=EX
			bpop=npop$pops
			N2=npop$N
			N2[EX]=b
			#step 3: simulate more generations
			npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=Time[j],r=r,K=K,n=n)
			}
		Exdat=cbind(Exdat,npop$sam)
		f=f+Fst(cbind(npop$pnames,npop$sam))
		}
	list(data=Exdat,Fst=f/L,extinctions=data.frame(island,Time))
	}


VI. Metapopulation Model- Extinction - Recolonization
# Function MPM is for a Metapopulation Model (MM) with extinction-recolonization
MPMer=function(r,K,L,m,mu,Time,n){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE	))
	f=matrix(0,ncol=12,nrow=12)
	events=length(Time)
	island=NA
	# first locus
	# Step 1: simulate for 10N generations to achieve equilibrium
	npop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)
	# Step 2: a random population goes extinct at each Time
	for(j in 1:events){
		EX=j
		island[j]=EX
		bpop=npop$pops
		bpop[EX,]=0
		N2=npop$N
		N2[EX]=1
		#Step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=Time[j],r=r,K=K,n=n)
		}
	Exdat=cbind(npop$pnames,npop$sam)
	f=f+Fst(cbind(npop$pnames,npop$sam))
	for(i in 2:L){
		# Step 1: simulate for 10N generations to achieve equilibrium
		npop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)
		# Step 2: a random population goes extinct at each Time
		for(j in 1:events){
			EX=j
			island[j]=EX
			bpop=npop$pops
			bpop[EX,]=0
			N2=npop$N
			N2[EX]=1
			#Step 3: simulate more generations
			npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=Time[j],r=r,K=K,n=n)
			}
		Exdat=cbind(Exdat,npop$sam)
		f=f+Fst(cbind(npop$pnames,npop$sam))
		}
	list(data=Exdat,Fst=f/L,extinctions=data.frame(island,Time))
	}

VII. Cornuet and Luikart Model
# function for 12 simultaneous bottlenecks and recovery to arbitrary size
Bottle12=function(r,K,L,m,mu,time,n,b,r2,K2){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		bpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 drops to b individuals
		N2=rep(b,12)
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r2,K=K2,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

Step 3) Migration Matrices

# No Migration Matrix
# simulate isolated populations	
m=0
x=c(
1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,
.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,
.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,
0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,
0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,
0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,
0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,
0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,
0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,
0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,
.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,
.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m)
m=matrix(x,nrow=12,ncol=12)

#Migration Matrix
#Adjust m for migration amount relative to original effective population size
#Example is for Nm = 1, Ne(o) = 1000
m=1/1000
x=c(
1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,
.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,
.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,
0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,
0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,
0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,
0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,
0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,
0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,
0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,
.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,
.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m)

m=matrix(x,nrow=12,ncol=12)


Step 4) Simulation File- input parameters and # of simulations
I. Focal Deme and Metapopulation Models
#Metasim runs the simulations for the following Population Models: FDMer, FDMbot, MPMer, or MPMbot
#Input is as follows:
#s=number of sims
#r=intrinsic rate of growth
#K=Effective population size (carrying capacity) for each deme
#L=number of loci drawn from each deme
#m=user defined migration matrix
#mu=mutation rate
#n=number of samples from each deme
#MM=the metapopulation model (either FDMer, FDMbot, MPMer, or MPMbot)

Metasim=function(s,r,K,L,m,mu,n,MM){
#Open blank files: for genotype file, Fst file, and event times file
write.table(NA,file="Genotypes.txt")
write.table(NA,file="Fst.txt")
write.table(NA,file="Times.txt")
for(i in 1:s){
	Times=rexp(10,1/500)
	Times=round(Times)
	Times=Times[order(Times)]
	Times=c(Times[1],diff(Times))
	Times=replace(Times,Times<1,2) #Time=0 gives 2 generations
	Sims=MM(r,K,L,m=m,mu,Time=Times,n)
	#Fill Data Tables
	write.table(Sims$data,file="Genotypes.txt",append=TRUE)
	write.table(Sims$Fst,file="Fst.txt",append=TRUE)
	write.table(Sims$extinctions,file="Times.txt",append=TRUE)
	}
}
#Example:
#Metasim(200,0.0002,1000,20,m,1/1000,20,MPMb)


II. Focal Deme Model with varying Ne
#unequalNe runs the simulations for the following Population Model: FPMN
#Input is as follows:
#s=number of sims
#r=intrinsic rate of growth
#Ne=effective populaton size of all demes
#Nb=effective population size for Deme 1
#K=Effective population size (carrying capacity) for each deme
#L=number of loci drawn from each deme
#m=user defined migration matrix
#mu=mutation rate
#n=number of samples from each deme
unequalNe=function(s,r,Ne,Nb,L,m,mu,time,n){
for(i in 1:s){
	B1000=FPMN(r,K=Ne,N=Nb,L,m,mu,time,n,b=Nb)
	write.table(B1000$data,file="Geontypes.txt",append=TRUE)
	write.table(B1000$Fst,file="Fst.txt",append=TRUE)
	}
}

#Example:
#unequalNe(1,.001,10,1,20,m,1/100,10,20)


III. Cornuet and Luikart Model

#For Cornuet and Luikart Model use the simulation file below
#Input is as follows:
#r=intrinsic rate of growth
#Ne=effective populaton size of all demes
#Nb=effective population size for Deme 1
#L=number of loci drawn from each deme
#m=user defined migration matrix
#mu=mutation rate
#n=number of samples from each deme
#MM=the metapopulation model (either FDMer, FDMbot, MPMer, or MPMbot)
CL=function(r,Ne,Nb,L,m,mu,time,n,alpha=Ne/Nb){
#input tau
tau=0.25
times <- rep(tau*Nb,8) # not clear for Cornuet and Luikart whether this should be Ne or Nb
write.table(data.frame(alpha,tau,times),file="CL_Genotypes.txt")
write.table(data.frame(alpha,tau,times),file="CL_Fst.txt")

for(i in 1:length(times)){
	B1000=Bottle12(r,K=Ne,L,m,mu,time=times[i],n,b=Nb,r2=r,K2=Ne)
	write.table(B1000$data,file="CLR100.8.txt",append=TRUE)
	write.table(B1000$Fst,file="CLR100.8.Fst.txt",append=TRUE)
	}
	}
#For population recovery post-bottleneck, add a value >0 for r
#Example: CL(.001,100,8,20,m,1/100,times,20,alpha)
