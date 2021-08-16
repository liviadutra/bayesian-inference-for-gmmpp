#include <full_conditionals_epidem.ox>

mcmc_func(const y, const S, const S_pred, const lambda_pred,
			const e, const k, const b2, psi, const psi_tamanho,
			const dist_inicial, const Q, X, const trunc_W,
			const prioris_psi, const par_priori_psi,
			const iter, const burnin){

decl time1;
time1 = timer(); 

decl E, matriz_cov, calib, stepgrid, grid, i, j, i1; 

E = range(1,e) - 1;

matriz_cov = new array[k];
for(j=0; j<k; ++j){
	matriz_cov[j] = (diag(psi[j]-1) .* 0.1) .^2;
}


calib = 2000;  

stepgrid = 1000; 
grid = range(S[0],S[1]+S_pred, (S[1]+S_pred) ./ stepgrid);

decl const_calib, aceitPsi, total_testesPsi;

const_calib = ones(1,k);

aceitPsi = zeros(1, k);	  
total_testesPsi = zeros(1, k);

decl Omega, B_matrix;
Omega = 5*diagonal(-Q)';			  
B_matrix = unit(e,e) + (1 ./ Omega) .* Q;	


decl W, X_aux, obs_bloco, psi_aux, psi_aux2;

W = cond_W(X, S, Q); 
W = W[0][] | (W[vecindex(W .>=trunc_W[0] .&& W .<=trunc_W[1])][]); 

X = cond_V(psi, y, W, B_matrix, S, b2);
X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);  
obs_bloco = countc(y,X[][1]);
psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi, &total_testesPsi,
						prioris_psi, par_priori_psi, b2)[0];							
psi_aux2 = <>;
for(j = 0; j < k; ++j){
	psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j]; 
	psi_aux2 ~= psi_aux[j];
}
		 			
decl const_calib_vec, probs_aceit_vec, prob_limite, i2;

const_calib_vec = zeros(1,k) | 2; 

decl iterCalib_MCMC;

iterCalib_MCMC = calib_MCMC(calib, X, S, Q, psi, y, Omega, B_matrix, k, e, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, prioris_psi, par_priori_psi, trunc_W, b2);	 
i2 = 1;

X = iterCalib_MCMC[0];	
psi = iterCalib_MCMC[1];	
psi_aux2 = iterCalib_MCMC[2][calib-1][];

while(sumr(const_calib_vec[sizer(const_calib_vec)-1][] .!=
			const_calib_vec[sizer(const_calib_vec)-2][]).!=0){	 
			
decl aux;
aux = cumulate(psi_tamanho' )';	
for(i1 = 0; i1 < k; ++i1){	
	matriz_cov[i1] = variance(iterCalib_MCMC[2][][(aux[i1] - psi_tamanho[i1]):
								(aux[i1] - 1)]) .* const_calib[i1];  
	matriz_cov[i1] = mpd(matriz_cov[i1], 10^(-5)); 
} 

aceitPsi = zeros(1, k);	  
total_testesPsi = zeros(1, k);		

iterCalib_MCMC = calib_MCMC(calib, X, S, Q, psi, y, Omega, B_matrix, k, e, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, prioris_psi, par_priori_psi, trunc_W, b2);
i2 += 1;
				
X = iterCalib_MCMC[0];	
psi = iterCalib_MCMC[1];	
psi_aux2 = iterCalib_MCMC[2][calib-1][]; 
				 
const_calib = calibracao(const_calib, aceitPsi, total_testesPsi);
const_calib_vec |= const_calib;	 
						  
}
  
decl psi_post, logposteriori, integral;

psi_post = constant(.NaN, burnin+iter, sumr(psi_tamanho));
logposteriori = constant(.NaN,burnin+iter,1);
integral = constant(.NaN, burnin+iter,1);

aceitPsi = zeros(1, k);	  
total_testesPsi = zeros(1, k);

for(i = 0; i < burnin; ++i){

	W = cond_W(X, S, Q);
	W = W[0][] | (W[vecindex(W .>=trunc_W[0] .&& W .<=trunc_W[1])][]);
	X = cond_V(psi, y, W, B_matrix, S, b2); 
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);  
	logposteriori[i] = cond_theta(k,e,X_aux,Q);	
	obs_bloco = countc(y,X[][1]);
	psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi, &total_testesPsi,
								prioris_psi, par_priori_psi, b2)[0];		
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j]; 
		psi_aux2 ~= psi_aux[j];	
	}	
											
	psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi, &total_testesPsi,
							prioris_psi, par_priori_psi, b2);
	logposteriori[i] +=	psi_aux[1];	
	psi_aux = psi_aux[0];
	
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j]; 
		psi_aux2 ~= psi_aux[j];	
	}	

	psi_post[i][] = psi_aux2;
	
	logposteriori[i] +=	logverossimilhanca(psi,y, obs_bloco, X, S, b2)[0];

	integral[i] = logverossimilhanca(psi,y, obs_bloco, X, S, b2)[1];						

}
	
decl X_iter, lambda_iter, int_cred, ajuste;

lambda_iter = constant(0, 1, columns(grid));
int_cred = constant(0,iter,(stepgrid+1));
	
aceitPsi = zeros(1, k);	  
total_testesPsi = zeros(1, k);

decl tempo_if, integral_previsao;
tempo_if = constant(0,iter,1);		  
integral_previsao = constant(0,iter,1);

for(i = 0; i < iter; ++i){

   	W = cond_W(X, S, Q);
	W = W[0][] | (W[vecindex(W .>=trunc_W[0] .&& W .<=trunc_W[1])][]);
	
	X = cond_V(psi, y, W, B_matrix, S, b2); 

	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);  
	logposteriori[burnin+i] = cond_theta(k,e,X_aux,Q);	
	obs_bloco = countc(y,X[][1]);
	psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi, &total_testesPsi,
							prioris_psi, par_priori_psi, b2)[0];		
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j]; 
		psi_aux2 ~= psi_aux[j];	
	}	
									
	psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi, &total_testesPsi,
							prioris_psi, par_priori_psi, b2);
	logposteriori[burnin+i] +=	psi_aux[1];
	psi_aux = psi_aux[0];
				
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j]; 
		psi_aux2 ~= psi_aux[j];	
	}	
	
	decl aux1,lambda;
	aux1 = sumc(X[][1] .<= grid)-1;										  
			
	lambda = func_g(psi, grid[vecindex(grid .< X[1][1])]', X[1][1], b2)[][0] |
			func_g(psi, grid[vecindex(grid .>= X[1][1])]', X[1][1], b2)[][1];	 
	X_iter = X[aux1][0]';	
	lambda_iter += lambda';	 	  
	int_cred[i][] = lambda'[range(0,stepgrid)];  	

	psi_post[burnin+i][] = psi_aux2;	 

	logposteriori[burnin+i] += logverossimilhanca(psi,y, obs_bloco, X, S, b2)[0];
	integral[burnin+i] = logverossimilhanca(psi,y, obs_bloco, X, S, b2)[1];
 

	if(fmod(i+1,10000)==0){println("time for ", burnin+i+1," iter: ",timespan(time1));}	
					
	decl previsao;

	previsao = predicao(psi,X[1][1],S, lambda_pred, b2);
							
	tempo_if[i][] = previsao[0];		  
	integral_previsao[i][] = previsao[1];
	
}

	println("accept probability Psi:", aceitPsi ./ total_testesPsi);

	savemat("lambda_mean.mat", lambda_iter ./ iter, 1);
	savemat("lambda_ci.mat", quantilec(int_cred,<0.025,0.975>), 1);
	savemat("psi.mat",psi_post,1);
	savemat("logposterior.mat",logposteriori,1);
	savemat("integral.mat",integral,1);
	savemat("timepredic.mat", tempo_if, 1);	
	savemat("integralpredic.mat", integral_previsao, 1);


}