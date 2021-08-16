#include <oxstd.h>
#include <oxprob.h>

#include <GMMPP_fixed_intercept.ox>


main(){
			

decl y, S, S_pred, e, k;
y = loadmat("data.mat");  // Data: n x 1. 				
S = <0,50>;	// observed time interval
S_pred = 0;	// prediction time

e = 3; // length of the ctmc state space	 
k = 3; // number of g functions


// Psi initial values
decl
a1 = <2.5>,					
b1 = <1.25>,
a2 = <30>,
b2 = <-1.25>,
a3 = <2.5>,
psi = {a1 ~ b1, a2 ~ b2, a3};

decl  dist_inicial, Q, X;
dist_inicial = <1,1,1> ./ 3; // initial distribution

Q = <-1/20, 1/20*1/2, 1/20*1/2; // initial Q-matrix
     1/20*1/2, -1/20, 1/20*1/2;
     1/20*1/2, 1/20*1/2, -1/20>;

X = <2,0>; // initial ctmc path
														
decl gamma, dir, gamma_psi, prioris_psi, par_priori_psi, psi_tamanho, g_nao_constante, g_constante; 
gamma = <>; // prior dist parameters for \theta. (diagonal): |K| x 2.
			// Or "= <>" (empty) if theta is fixed.
dir = <>; // prior dist parameters for \theta.. : |K| x |K|
		  // Or "= <>" (empty) if theta is fixed.

gamma_psi = {<1,1>};	// prior dist parameters for constant \psi: array[i] = 1 x 2.
 
prioris_psi = {"uniform","uniform","uniform","uniform"}; // prior distributions for psi's nonconstant
														 // Options: "uniform" or "normal"
par_priori_psi = {<>,<>,<>,<>};	// prior dist parameters. For "uniform": <>;
								//						  For "normal": <mean,variance>. 
		
psi_tamanho = <2, 2, 1>; // vector: i-th column is the number of parameters in the functional form  g_i
g_nao_constante = <0,1>; // which position of "psi" are nonconstant functional forms
g_constante = <2>; // which position of "psi" are constant functional forms


decl iter, burnin, iterWV;
iter = 90000;
burnin = 10000;
iterWV = 5;

decl limite_trunc_exp;
limite_trunc_exp = {+.Inf, - psi[1][0] ./ psi[1][1], +.Inf};
					// for prediction: this variable defines the time limit
	                // for decreasing functional forms. Negative IF cannot be generated.

mcmc_func(y, S, S_pred, e, k, psi, dist_inicial, Q, X, gamma, dir, gamma_psi, prioris_psi,
			par_priori_psi, psi_tamanho, g_nao_constante, g_constante, iter, burnin, iterWV, limite_trunc_exp);

}																	
