#include <oxstd.h>
#include <oxprob.h>

#include <GMMPP_var_intercept.ox>


main(){

decl y, S, S_pred, e, k;
y = loadmat("data.mat"); // Data: n x 1. 	
S = <0,112>;  	// observed time interval
S_pred = 0;	// prediction time

e = 3;	// length of the ctmc state space	 
k = 3;	// number of g functions

// Psi initial values
decl
a1 = <3>,
b1 = <-0.1>,
a2 = <0.5>,
b2 = <0.1>,
psi = {b1, b2};

decl dist_inicial, Q, X, R;
dist_inicial = <1,1,1> ./ 3;	// initial distribution					

// initial Q-matrix
Q = <-1/60, 1/120, 1/120;
	  1/120,-1/60, 1/120;
	  1/120, 1/120, -1/60>;

// initial ctmc path	  
X = <2,0;
	 0,30;
	 2,50>;	
	 
R = <3.5,3.5,1>'; // initial R values.
							// = 0 for which the functional form does not have variable intercept

// defining prior for R|V=k, for each g_k
decl  classe_priori_r, valores_r, probs_r, hiperpar_r;

classe_priori_r = {"R_discreta","R_discreta","R_uniforme"};		// Options: "R_discreta", "R_uniforme" or "R_gamma"
valores_r = {range(0.1,4.7,0.2)', range(0.1,4.7,0.2)', 0}; // if classe_priori_r[i] = R_uniforme or R_gamma, then, valores_r[i] = <0>;
probs_r = { ones(1,24)./24,  ones(1,24)./24, ones(1,24)./24};	    // if classe_priori_r[i] = R_uniforme or R_gamma, then, probs_r[i] = <1>;

hiperpar_r = {<>,<>,<>};		 // if classe_priori_r[i] = R_gamma, then hiperpar_r[i] = <alpha, beta>
							 // if classe_priori_r[i] = R_uniforme or R_discreta, then, hiperpar_r[i] = <>; 


decl gamma, dir, gamma_psi, prioris_psi, par_priori_psi, psi_tamanho, g_nao_constante, g_constante; 
gamma = <>;	// prior dist parameters for \theta. (diagonal): |K| x 2.
						// Or "= <>" (empty) if theta is fixed.
dir = constant(1, k, k);  // prior dist parameters for \theta.. : |K| x |K|
		  				  // Or "= <>" (empty) if theta is fixed.

gamma_psi = {<>};	 	// prior dist parameters for constant \psi: array[i] = 1 x 2.  

prioris_psi={"normal","normal"}; // prior distributions for psi's nonconstant
								   // Options: "uniform" or "normal"
par_priori_psi = {<-0.1,0.2^2>,<0.1,0.2^2>};	// prior dist parameters. For "uniform": <>;
							// For "normal": <mean,variance>. 
		


psi_tamanho = <1, 1>; // vector: i-th column is the number of parameters in the functional form  g_i
g_nao_constante = <0,1>; // which position of "psi" are nonconstant functional forms
g_constante = <>; // which position of "psi" are constant functional forms




decl iter, iterWVR, burnin;

iter = 90000;
iterWVR = 1;
burnin = 10000;


mcmc_func(y, S, S_pred, e, k, psi, dist_inicial, Q, X, R,
			gamma, dir, gamma_psi, valores_r, probs_r, classe_priori_r, hiperpar_r,
			prioris_psi, par_priori_psi, psi_tamanho,
			 g_nao_constante, g_constante, iter, burnin, iterWVR);
 				
}
