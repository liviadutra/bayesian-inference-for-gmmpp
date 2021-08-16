#include <full_conditionals.ox>

mcmc_func(const y, const S, const S_pred, const e, const k,
			psi, const dist_inicial, const Q, X,
			const gamma, const dir, const gamma_psi, const prioris_psi,
			const par_priori_psi, const psi_tamanho, const g_nao_constante,
			const g_constante, const iter, const burnin, const iterWV, 
            const limite_trunc_exp){

decl time1;
time1 = timer(); 

decl E, matriz_cov, calib, stepgrid, grid, i, j, i1; 

E = range(1,e) - 1;

matriz_cov = new array[columns(g_nao_constante)];
for(j=0; j<columns(g_nao_constante); ++j){
	matriz_cov[j] = (diag(psi[g_nao_constante[j]]) .* 0.1) .^2;
}				

calib = 1000;  // numero de iterações a cada ciclo para calibração

stepgrid = 1000; // number of points to plot the estimated IF - mean and CI
grid = range(S[0],S[1], S[1] ./ stepgrid);


decl const_calib, num_intervalos, aceitCTMC, total_testesCTMC, aceitPsi, total_testesPsi;

const_calib = ones(1,columns(g_nao_constante));
num_intervalos = 1;

aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));

decl Omega, B_matrix;
Omega = 2*diagonal(-Q)';			  
B_matrix = unit(e,e) + (1 ./ Omega) .* Q;	


// primeira iteração:
decl W, X_aux, Q_matrix, intervalos_tempos, obs_bloco, psi_aux, psi_aux2;

for(j = 0; j < iterWV; ++j){
	W = cond_W(X, S, Q); 		  
	X = cond_V(psi, num_intervalos, y, W, X, S, B_matrix, dist_inicial, &aceitCTMC, &total_testesCTMC);
}

X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);
Q_matrix = cond_theta(k, e, X_aux, Q, gamma, dir)[0];	  								 
Omega = 2*diagonal(-Q_matrix)';			  	
B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];	

obs_bloco = countc(y,X[][1]);			
psi_aux = cond_Psi_MH(e, X ~ X_aux[][1], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi,
							prioris_psi, par_priori_psi)[0];							
psi_aux2 = <>;
for(j = 0; j < k; ++j){
	psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];
	psi_aux2 ~= psi_aux[j];
}							  

				   
decl num_intervalos_vec, const_calib_vec, probs_aceit_vec, prob_limite, i2;

num_intervalos_vec = 0 | 10; 
const_calib_vec = zeros(1,columns(g_nao_constante)) | 2;

probs_aceit_vec = new array[0];	 
prob_limite = 0.25;				


decl iterCalib_MCMC;

iterCalib_MCMC = calib_MCMC(calib, iterWV, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2, prioris_psi, par_priori_psi);
i2 = 1;

X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][];

while(num_intervalos_vec[sizer(num_intervalos_vec)-1] .!=
			num_intervalos_vec[sizer(num_intervalos_vec)-2] ||
		sumr(const_calib_vec[sizer(const_calib_vec)-1][] .!=
			const_calib_vec[sizer(const_calib_vec)-2][]).!=0){

decl aux;
aux = cumulate(psi_tamanho' )';	
for(i1 = 0; i1 < columns(g_nao_constante); ++i1){
	matriz_cov[i1] = variance(iterCalib_MCMC[3][][aux[g_nao_constante[i1]] - psi_tamanho[g_nao_constante[i1]]:
								aux[g_nao_constante[i1]] - 1]) .* const_calib[i1];
	matriz_cov[i1] = mpd(matriz_cov[i1], 10^(-5));   
}

aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));		

iterCalib_MCMC = calib_MCMC(calib, iterWV, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2, prioris_psi, par_priori_psi);
i2 += 1;
				
X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][]; 

probs_aceit_vec |= {aceitCTMC ./ total_testesCTMC};		 

if(quantiler(aceitCTMC ./ total_testesCTMC) <= prob_limite){   
	num_intervalos += 1;	   
}
	
aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));		

iterCalib_MCMC = calib_MCMC(calib, iterWV, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2, prioris_psi, par_priori_psi);
i2 += 1;
				
X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][]; 
				 
const_calib = calibracao(const_calib, aceitPsi, total_testesPsi); 													  
num_intervalos_vec |= num_intervalos;
const_calib_vec |= const_calib;
						  
}											
	
// BURNIN
		  
decl theta_post, psi_post, logposteriori, integral, tamanho_X, tempo_salto1;

theta_post = constant(.NaN, burnin+iter, e*e);
psi_post = constant(.NaN, burnin+iter, sumr(psi_tamanho));
logposteriori = constant(.NaN,burnin+iter,1);
integral = constant(.NaN, burnin+iter,1);
tamanho_X = constant(0,burnin+iter,1);
tempo_salto1 = constant(0,burnin+iter,1);

aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));

for(i = 0; i < burnin; ++i){

	for(j = 0; j < iterWV; ++j){
		W = cond_W(X, S, Q_matrix); 	
		X = cond_V(psi, num_intervalos, y, W, X, S, B_matrix, dist_inicial, &aceitCTMC, &total_testesCTMC);
	}

	tamanho_X[i] = rows(X);
	tempo_salto1[i] = X[1][1];
	
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);
	Q_matrix = cond_theta(k, e, X_aux, Q_matrix, gamma, dir);
	logposteriori[i] = Q_matrix[1];	 
	Q_matrix = Q_matrix[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];
	obs_bloco = countc(y,X[][1]);
	
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi,
							prioris_psi, par_priori_psi);
	logposteriori[i] +=	psi_aux[1];	
	psi_aux = psi_aux[0];

	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];
		psi_aux2 ~= psi_aux[j];
	}

	theta_post[i][] = vecr(Q_matrix)';
	psi_post[i][] = psi_aux2;
	
	logposteriori[i] +=	logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, range (0, rows(X)-1))[0]+
							log(dist_inicial[X_aux[0][0]]);
	integral[i] = sumc(logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, range (0, rows(X)-1))[1]);						

}	
	
// INICIO DA CADEIA

decl X_iter, lambda_iter, int_cred;

lambda_iter = constant(0, 1, columns(grid));
int_cred = constant(0,iter,stepgrid+1);
	
aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));

decl PP_previsao, integral_previsao;
PP_previsao = constant(0,iter,1);
integral_previsao = constant(0,iter,1);
	 
for(i = 0; i < iter; ++i){

	for(j = 0; j < iterWV; ++j){
		W = cond_W(X, S, Q_matrix); 
		X = cond_V(psi, num_intervalos, y, W, X, S, B_matrix, dist_inicial, &aceitCTMC, &total_testesCTMC);
	}

	tamanho_X[burnin+i] = rows(X);
	tempo_salto1[burnin+i] = X[1][1];
	
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);  	
	Q_matrix =  cond_theta(k, e, X_aux, Q_matrix, gamma, dir);
	logposteriori[burnin+i] = Q_matrix[1];
	Q_matrix = Q_matrix[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];
	obs_bloco = countc(y,X[][1]);
	
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi, prioris_psi, par_priori_psi);
	logposteriori[burnin+i] += psi_aux[1];
	psi_aux = psi_aux[0];
	
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];
		psi_aux2 ~= psi_aux[j];
	}
	
	// guardando resultados
	decl aux1, aux2, aux3, lambda;
	aux1 = sumc(X[][1] .<= grid)-1;										  
	aux3 = <>;
	for(j = 0;j<e;++j){
		aux3 ~= X[aux1][1];
	}		
	lambda = func_g(psi, grid',aux3);	 
	X_iter = X[aux1][0]';	
	aux2 = selectrc(lambda,range(0,columns(grid)-1, 1)', X_iter);  
	lambda_iter += aux2;	 
	int_cred[i][] = aux2[range(0,stepgrid)];
	
	theta_post[burnin+i][] = vecr(Q_matrix)';
	psi_post[burnin+i][] = psi_aux2;	 

	logposteriori[burnin+i] += logverossimilhanca(e,psi,y,obs_bloco,X,intervalos_tempos,range(0,rows(X)-1))[0]+
							   log(dist_inicial[X_aux[0][0]]);
	integral[burnin+i] = sumc(logverossimilhanca(e,psi,y,obs_bloco,X,intervalos_tempos,range(0,rows(X)-1))[1]);
							   
	if(fmod(i+1,10000)==0){println("time for ", burnin+i+1," iter: ",timespan(time1));}

	if(S_pred > 0){
	decl previsao;

	previsao = gerar_processo(S, S_pred, E, Q_matrix, psi, e, 
							range(S[1],S[1]+S_pred,S_pred./1000),X[rows(X)-1][],
							limite_trunc_exp);
	PP_previsao[i][] = previsao[1];
	integral_previsao[i][] = previsao[2];
	}
	
}

	println("accept probability (V_0, V):", aceitCTMC ./ total_testesCTMC);
	println("accept probability Psi:", aceitPsi ./ total_testesPsi);

	savemat("lambda_mean.mat", lambda_iter ./ iter, 1);
	savemat("lambda_ci.mat", quantilec(int_cred,<0.025,0.975>), 1);
	savemat("psi.mat",psi_post,1);
	savemat("logposterior.mat",logposteriori,1);
	savemat("integral.mat",integral,1);
	savemat("CTMCjumps.mat",tamanho_X,1);
	savemat("firstjump.mat",tempo_salto1,1);

	if(rows(gamma) .>= 1){
	savemat("theta.mat",theta_post,1);
	}

	if(S_pred > 0){	 
	savemat("Npredic.mat", PP_previsao, 1);	
	savemat("integralpredic.mat", integral_previsao, 1);
	}

}
