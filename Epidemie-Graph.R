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
	Markov=c(0,0.75,0.24,0.01) 					#Chaine de Markov pour les individus infectés (état 2)
	E=matrix(nrow=n,ncol=500)				    #Matrice des états de chaque individu, à chaque instant t.
	E[,1]=1 										#On initialise tous les individus à l'état 1
	E[,2:500]=0 								#On remplit les autres colones avec des 0 (cf ligne41)
	IMN=sample(n,floor(n/10))				    #Un dizième de la pop est Immunisée naturelement
	E[IMN,1:500]=3								#On applique l'état 3 à cette part de la population jusqu'à la fin.
	P0=sample(1:n,1) 							#On séléctionne le patient 0
	E[P0,1]=2									#On ajoute P0 à l'état infecté au temps 1

	plot(X[which(E[,1]==1)],Y[which(E[,1]==1)], pch = 20,xlab='X',ylab='Y',xlim=c(-4.5,4.5),ylim=c(-4.5,4.5), main= 'temps : t=0') 
	points(X[which(E[,1]==3)],Y[which(E[,1]==3)], pch = 20, col ='green') 
	points(X[which(E[,1]==2)],Y[which(E[,1]==2)], pch=20,col='red') 	
	T=1 
	while (2%in%E[,T]){ #on déroule la chaine tant qu'il y a encore des malades.
			
		SA=c(which(E[,T]==1)) 								#On isole les individus à l'état 1 au temps T dans un vecteur.
		IN=c(which(E[,T]==2)) 
		POS=matrix(c(X,Y),nrow=n) 							#Matrice des coordonnées
		DISTALL=matrix(as.matrix(dist(POS)),nrow=n) 			#Matrice des distances 
		DIST=as.matrix(DISTALL[SA,IN])						#Distances entre infectés et Sains
		PROCHE=which(apply(DIST,1,function(x) any(x<0.15))) #Retourne les indices des lignes de la matrice DIST 
		E[SA,T+1]=1
		E[SA[PROCHE],T+1]=2                                 #Applique l'état 2(infecté) aux individus sains correspondant aux indices de PROCHE /!\C'est ici que ça merde/!\
					
	
		E[IN,T+1]=sapply(E[IN,T+1],function(x) x+sample(1:4,1,prob=Markov)) #Chaine de markov pour les individus infectés 
		E[which(E[,T+1]==3),(T+2):500]=3  #On applique l'état 3 jusqu'à la fin pour les individus qui viennennt d'entrer dans cet état.
		E[which(E[,T+1]==4),(T+2):500]=4  #On applique l'état 4 jusqu'à la fin pour les individus qui viennennt d'entrer dans cet état.
	
		X=sapply(X,function(x) x+runif(1,-0.15,0.15)) #Déplacement des individus selon l'axe X
		Y=sapply(Y,function(x) x+runif(1,-0.15,0.15)) #Déplacement des individus selon l'axe Y
			
		plot.new()
		plot(X[which(E[,T+1]==1)],Y[which(E[,T+1]==1)], pch = 20,xlab='X',ylab='Y', xlim=c(-5,5),ylim=c(-5,5), main= T)
		points(X[which(E[,T+1]==3)],Y[which(E[,T+1]==3)], pch = 20, col ='green') 
		points(X[which(E[,T+1]==2)],Y[which(E[,T+1]==2)], pch=20,col='red')
		points(X[which(E[,T+1]==4)],Y[which(E[,T+1]==4)], pch=20,col='blue')
		T=T+1	
		next
	}
	return(E[,1:T])
}

