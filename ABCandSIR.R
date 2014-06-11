# Modèle SIR
# Pop fermée & Susceptible -> Infectious -> Recovered
# S_t,I_t redonnent R_t
# Tk= T_k-1 + exp

#modèle SIR classique

SIR=function(N,lambda,gamma){
  T=cbind( cumsum(rexp(N,lambda)) , rep(1,N) ) ;    #dates d'infection
  D=cbind( T[,1]+rexp(N,gamma) , rep(-1,N) );       #dates de rémission

  
  x=rbind(T,D);
  I=x[order(x[,1]),]
  
  plot(I[,1],cumsum(I[,2]),type="s",col="blue",main="SIR",ylim=c(-10, N),ylab="",xlab="time"); #infectious
  lines(T[,1],rep(N,N)-1:N,type="s", col="red"); #saine
  DD=D[order(D[,1]),]
  lines(DD[,1],1:N,type="s",col="green") #recovered
	legend(5,80, c("S","I","R"),
	lty=c(1,1,1),lwd=c(2.5,2.5),col=c("red","blue","green"));
}

SIR(100,50,1);

#ajour du contact tracing

SIRCT=function(N,S0,I0,lambda1,lambda2,lambda3,mu,plot){
	T=rexp(N,1);
	C=rep(0,N);
	S=rep(0,N);
	I=rep(0,N);
	R1=rep(0,N);
	R2=rep(0,N);

	S[1]=S0;
	I[1]=I0;
	C[1]= lambda1*S[1]*I[1]+(mu+lambda2)*I[1];
	
	for (k in 2:N){
		U=runif(1,0,C[k-1]);
		S[k] = S[k-1]- 1*(U<lambda1* S[k-1] * I[k-1]);
		I[k] = I[k-1] - 1*(U< lambda3* I[k-1] )*(U>lambda1*S[k-1]*I[k-1]) + 1*(U<lambda1*S[k-1]*I[k-1]); # rajouter le coeff dans la somme
		R1[k] = R1[k-1]+1*(U>lambda1*S[k-1]*I[k-1]+mu*I[k-1])*(U<lambda1*S[k-1]*I[k-1]+(mu+lambda2)*I[k-1]);
		R2[k] = R2[k-1]+1*(U>lambda1*S[k-1]*I[k-1]+(mu+lambda2)*I[k-1])*(U<lambda1*S[k-1]*I[k-1]+(mu+lambda2)*I[k-1]+lambda3*I[k-1]);
		C[k] = lambda1 *S[k-1]*I[k-1] + (mu+lambda2)*I[k-1] + lambda3*I[k-1]*(R1[k-1]+R2[k-1]);		
		}  
	
	TT= ((C+0.001)^(-1))*T;

	if (plot==TRUE){
	plot(cumsum(TT),S,type="s",col="blue",main="SIR with contact tracing",ylim=c(-10, S0+I0),ylab="Process",xlab="Time") #saine
  	lines(cumsum(TT),I, type="s", col="red") #infectious
	lines(cumsum(TT), R1, type="s", col="green");
	lines(cumsum(TT), R2, type="s", col="yellow");

	legend(2,205, c("S","I","R1","R2"),
	lty=c(1,1,1,1),lwd=c(2.5,2.5),col=c("blue","red","green","yellow")); }

	return(cbind(cumsum(TT),R1,R2));
}

SIRCT(10000,200,1,0.5,0.2,0.3,0.5,TRUE);

#méthode abc

library(abc);
Robs=SIRCT(10000,200,1,0.5,0.2,0.3,0.5,FALSE);

N=10;

#simulations des paramètres

mu=runif(N,-6,-4);
lambda1=runif(N,-9,-6);
lambda2=runif(N,-4,3);
lambda3=runif(N,-8,2);

TTsim=matrix(0,10000,N);
R1sim=matrix(0,10000,N);
R2sim=matrix(0,10000,N);

for (j in 1:N){
	S=SIRCT(10000,200,1,exp(lambda1),exp(lambda2),exp(lambda3),exp(mu[j]),FALSE);
	TTsim[,j]=S[,1];
	R1sim[,j]=S[,2];
	R2sim[,j]=S[,3];
}

Tobs=cbind( Robs[,1] , rep(1,N) ) ;      
Tsim=cbind( TTsim[,1] , rep(-1,N) );      
  
x=rbind(Tobs,Tsim);
I=x[order(x[,1]),];











