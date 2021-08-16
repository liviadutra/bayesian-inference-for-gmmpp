#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.oxh>
#include <GMMPP_epidem.ox>

 
main(){
								   
decl y, S, S_pred, lambda_pred, e, k;
y = loadmat("data.mat");  // Data: n x 1.				
S = <0,95>;	   // observed time interval
S_pred = 65;	// prediction time
lambda_pred = 10; // the target value that the IF reaches after S

e = 2; // length of the ctmc state space
k = 2; // number of g functions


decl
b2 = 0;	// fixed value

// Psi initial values
decl
a1 = <360>,					
b1 = <0>,
c1 = <20>,
d1 = <-4>,			  
c2 = <8>,
d2 = <2.3>,
psi = {a1 ~ b1 ~ c1 ~ d1, c2 ~ d2},
psi_tamanho = <4,2>; 

decl dist_inicial, Q, X, trunc_W;
dist_inicial = <1,0> ./ 1;	// initial distribution					

Q = <-1/50, 1/50;	 // initial Q-matrix
     1/40, -1/40>;	

X = <0,0;	 // initial ctmc path
	 1,45>;
	 
trunc_W = <30,75>; // truncation of the change time


decl prioris_psi, par_priori_psi;
prioris_psi = {"uniform","uniform","uniform","uniform","uniform","uniform"};	 
								// prior distributions for psi's nonconstant
								// Options: "uniform" or "normal"
								
par_priori_psi = {<>,<>,<>,<>,<>,<>};
					// prior dist parameters. For "uniform": <>;
					//						  For "normal": <mean,variance>. 
				   


decl iter, burnin;
iter = 450000;
burnin = 50000; 


mcmc_func(y, S, S_pred, lambda_pred, e, k, b2, psi, psi_tamanho,
			dist_inicial, Q, X, trunc_W, prioris_psi, par_priori_psi,
			iter, burnin);
			
}																	
