epidemie=function(n){ 
	f=function (x){	  #fonction de densité de la loi mélange de deux lois normales.  
		coef=runif(1) #On tire le "poids" des deux lois normales aléatoirement
		coef*(1/sqrt(2*pi*0.5))*exp(-(x+2)^2) + (1-coef)*(1/sqrt(2*pi*0.5))*exp(-(x-2)^2)}
	F=function (x){	
		coef=runif(1)
		coef*(1/sqrt(2*pi*0.5))*exp(-(x+2)^2) + (1-coef)*(1/sqrt(2*pi*0.5))*exp(-(x-2)^2)}
	x=seq(-5,5,by=0.0001)
	X=c() #vecteur abscisse
	Y=c() #vecteur ordonnée
	X=sample(x,n,rep=TRUE,prob=f(x))			#On simule n abscisses selon notre loi de densité
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
	Z1=apply(DISTALL,1,function(x) sum(x<0.4))				#On compte le nombre d'individus dans la zone 1 de chaque individu
	Z2=apply(DISTALL,1,function(x) sum(x<2))
	Z3=apply(DISTALL,1,function(x) sum(x>=2))
		
	T=1 
	while (2%in%E[,T]){ 									#on déroule la chaine tant qu'il y a encore des malades.
			
		SA=c(which(E[,T]==1)) 								#On isole les individus à l'état 1 au temps T dans un vecteur.
		IN=c(which(E[,T]==2)) 
		
		DIST=as.matrix(DISTALL[SA,IN])						#Distances entre infectés et Sains
		
		ZONE=matrix(nrow=length(SA),ncol=6)
		ZONE[,1]=apply(DIST,1,function(x) sum(x<0.4))   #Nb d'individus infectés zone 1
		ZONE[,2]=Z1[SA] 									#Nb d'individus zone 1
		ZONE[,3]=apply(DIST,1,function(x) sum(x<2))		#Nb d'individus infectés zone 2
		ZONE[,4]=Z2[SA]	 							    #Nb d'individus zone 2
		ZONE[,5]=apply(DIST,1,function(x) sum(x>=2))	
		ZONE[,6]=Z3[SA]
		
		p=0.8*ZONE[,1]/(ZONE[,2]+1)+0.15*ZONE[,3]/(ZONE[,4]+1)+0.05*ZONE[,5]/(ZONE[,6]+1) #Vecteur des probabilités de devenir infectés pour chaque individu sain
		E[SA,T+1]=sapply(p,function(x) sample(1:2,1,prob=c(1-x,x))) #Chaine de Markov pour les individus sains
		
		E[IN,T+1]=sapply(E[IN,T+1],function(x) x+sample(1:4,1,prob=Markov)) #Chaine de markov pour les individus infectés 
		E[which(E[,T+1]==3),T+2]=3  #On applique l'état 3 jusqu'à la fin pour les individus qui viennent d'entrer dans cet état.
		E[which(E[,T+1]==4),T+2]=4  #On applique l'état 4 jusqu'à la fin pour les individus qui viennent d'entrer dans cet état.
					
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
	ICt1=c()
	ICt2=c()
	ICm1=c()
	ICm2=c()

	for (i in 1:n){
		Ep=epidemie(t)
		X[i]=Ep[1]    #Vecteur des temps de disparitions de l'épidémie
		Y[i]=(Ep[2]+0.25*(Ep[3]-10))/2	  #Réduction de variance sur le pourcentage de décès.
		M[i]=mean(X[1:i])
		m[i]=mean(Y[1:i])
		S[i]=sd(X[1:i])
		s[i]=sd(Y[1:i])
		ICt1[i]=M[i]-1.96*S[i]/sqrt(i) 
		ICt2[i]=M[i]+1.96*S[i]/sqrt(i)
		ICm1[i]=m[i]-1.96*s[i]/sqrt(i)
		ICm2[i]=m[i]+1.96*s[i]/sqrt(i)
	}	

	ICt=ICt2[n]-ICt1[n] #Taille de l'intervalle de confiance pour l'espérance du temps de disparition.
	ICm=ICm2[n]-ICm1[n]
	par(mfrow=c(1,2))
	matplot(1:n,matrix(c(M[1:n],ICt1[1:n],ICt2[1:n]),ncol=3),ylim=c(0,60),col=c("black","red","red"),type='l',main="Convergence du temps de disparition de l'épidémie",xlab="Nombre de simulations",ylab="Temps de disparition moyen")
	matplot(1:n,matrix(c(m[1:n],ICm1[1:n],ICm2[1:n]),ncol=3),ylim=c(0,25),col=c("black","red","red"),type='l',main="Convergence du pourcentage de décès",xlab="Nombre de simulations",ylab="Pourcentage de décès moyen")
	#hist(X,main="Histogramme du temps de disparition",xlab="Temps de disparition",ylab="Fréquence")
	#hist(Y,main ="Histogramme du pourcentage de décès",xlab="Pourcentage de décès",ylab="Fréquence")
	fin=Sys.time()
	temps=fin-deb
	return(c(mean(X),M[n],m[n],ICt,ICt1[n],ICt2[n],var(X),ICm,ICm1[n],ICm2[n],var(Y),temps))
}