epidemie=function(n){ 
	f=function (x){	
		coef=runif(1)
		coef*(1/sqrt(2*pi*0.5))*exp(-(x+2)^2) + (1-coef)*(1/sqrt(2*pi*0.5))*exp(-(x-2)^2)}
	F=function (x){	
		coef=runif(1)
		coef*(1/sqrt(2*pi*0.5))*exp(-(x+2)^2) + (1-coef)*(1/sqrt(2*pi*0.5))*exp(-(x-2)^2)}
	x=seq(-5,5,by=0.0001)
	X=c() #vecteur abscisse
	Y=c() #vecteur ordonnée
	X=sample(x,n,rep=TRUE,prob=f(x))
	Y=sample(x,n,rep=TRUE,prob=F(x))
	Markov=c(0,0.75,0.2,0.05) 					#Chaine de Markov pour les individus infectés (état 2)
	E=matrix(nrow=n,ncol=500)				    #Matrice des états de chaque individu, à chaque instant t.
	E[,1]=1 										#On initialise tous les individus à l'état 1
	E[,2:500]=0 								#On remplit les autres colones avec des 0 (cf ligne41)
	IMN=sample(n,floor(n/10))				    #Un dizième de la pop est Immunisée naturelement
	E[IMN,1:500]=3								#On applique l'état 3 à cette part de la population jusqu'à la fin.
	P0=sample(1:n,1) 							#On séléctionne le patient 0
	E[P0,1]=2									#On ajoute P0 à l'état infecté au temps 1
	POS=matrix(c(X,Y),nrow=n) 								#Matrice des coordonnées
	DISTALL=matrix(as.matrix(dist(POS)),nrow=n) 			#Matrice des distances 
	Z1=apply(DISTALL,1,function(x) sum(x<0.4))				
	Z2=apply(DISTALL,1,function(x) sum(x<2))
	Z3=apply(DISTALL,1,function(x) sum(x>=2))
		
	T=1 
	while (2%in%E[,T]){ #on déroule la chaine tant qu'il y a encore des malades.
			
		SA=c(which(E[,T]==1)) 								#On isole les individus à l'état 1 au temps T dans un vecteur.
		IN=c(which(E[,T]==2)) 
		
		DIST=as.matrix(DISTALL[SA,IN])						#Distances entre infectés et Sains
		
		ZONE=matrix(nrow=length(SA),ncol=6)
		ZONE[,1]=apply(DIST,1,function(x) sum(x<0.4))          #Nb d'individus infectés zone 1
		ZONE[,2]=Z1[SA] 		#Nb d'individus zone 1
		ZONE[,3]=apply(DIST,1,function(x) sum(x<2))			   #Nb d'individus infectés zone 2
		ZONE[,4]=Z2[SA]	    #Nb d'individus zone 2
		ZONE[,5]=apply(DIST,1,function(x) sum(x>=2))	
		ZONE[,6]=Z3[SA]
		
		p=0.8*ZONE[,1]/(ZONE[,2]+1)+0.15*ZONE[,3]/(ZONE[,4]+1)+0.05*ZONE[,5]/(ZONE[,6]+1)
		E[SA,T+1]=sapply(p,function(x) sample(1:2,1,prob=c(1-x,x)))
		
		E[IN,T+1]=sapply(E[IN,T+1],function(x) x+sample(1:4,1,prob=Markov)) #Chaine de markov pour les individus infectés 
		E[which(E[,T+1]==3),T+2]=3  #On applique l'état 3 jusqu'à la fin pour les individus qui viennennt d'entrer dans cet état.
		E[which(E[,T+1]==4),T+2]=4  #On applique l'état 4 jusqu'à la fin pour les individus qui viennennt d'entrer dans cet état.
					
		T=T+1	
		next
	}
	return(c(T,100*length(which(E[,T]==4))/n,100*length(which(E[,T]==3))/n))
}

MC=function(n,t){
	deb=Sys.time()
	Ep=c()
	X=c()
	M=c()
	Y=c()
	m=c()
	S=c()
	s=c()
	ICm1=c()
	ICm2=c()
	ICv1=c()
	ICv2=c()
	for (i in 1:n){
		Ep=epidemie(t)
		X[i]=Ep[2]
		Y[i]=(Ep[2]+0.25*(Ep[3]-10))/2
		M[i]=mean(X[1:i])
		m[i]=mean(Y[1:i])
		S[i]=sd(X[1:i])
		s[i]=sd(Y[1:i])
		ICm1[i]=M[i]-1.96*S[i]/sqrt(i)
		ICm2[i]=M[i]+1.96*S[i]/sqrt(i)
		ICv1[i]=m[i]-1.96*s[i]/sqrt(i)
		ICv2[i]=m[i]+1.96*s[i]/sqrt(i)
	}	
	ICv=ICv2[n]-ICv1[n]
	ICm=ICm2[n]-ICm1[n]
	par(mfrow=c(2,2))
	matplot(1:n,matrix(c(M[1:n],ICm1[1:n],ICm2[1:n]),ncol=3),ylim=c(0,25),col=c("black","red","red"),type='l',main="Convergence du pourcentage de décès",xlab="Nombre de simulations",ylab="Pourcentage de décès moyen")
	matplot(1:n,matrix(c(m[1:n],ICv1[1:n],ICv2[1:n]),ncol=3),ylim=c(0,25),col=c("black","red","red"),type='l',main="Convergence de la nouvelle variable",xlab="Nombre de simulations",ylab="Moyenne de la nouvelle ")
	hist(X,main ="Histogramme du pourcentage de décès",xlab="Pourcentage de décès",ylab="Fréquence")
	hist(Y,main="Histogramme de la nouvelle variable",xlab="Nouvelle variable",ylab="Fréquence")
	fin=Sys.time()
	temps=fin-deb
	return(c(M[n],m[n],ICm1[n],ICm2[n],ICv1[n],ICv2[n],var(X),var(Y),temps))
}
