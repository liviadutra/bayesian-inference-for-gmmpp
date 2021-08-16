#include <full_conditionals_var.ox>

mcmc_func(const y, const S, const S_pred, const e, const k,
			psi, const dist_inicial, const Q, X, R,
			const gamma, const dir, const gamma_psi,
			const valores_r, const probs_r, const classe_priori_r, const hiperpar_r,
			const prioris_psi, const par_priori_psi, const psi_tamanho,
			const g_nao_constante, const g_constante, const iter, const burnin, const iterWVR){

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

for(j = 0; j < iterWVR; ++j){
	W = cond_W(X, S, Q);
	intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];  // matriz em que cada linha são os intervalos entre saltos da ctmc(virtuais ou não)  [w_i, w_{i+1})						 
	obs_bloco = dropr(countc(y,W), rows(W));	  // conta o número de observações em cada intervalo [w_i, w_{i+1}). 
											  // primeira linha obrigatoriamente é zero.
	X = cond_VR(psi,num_intervalos,y,W, X ~ R, S,
		B_matrix, dist_inicial, &aceitCTMC, &total_testesCTMC,
		intervalos_tempos, obs_bloco,
		k, valores_r, probs_r, classe_priori_r, hiperpar_r); 
}
		
X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);  	// X_aux é a trajetoria de X com tempo de permanencia em cada estado,
																// ao inves do tempo de salto.
Q_matrix = cond_theta(k, e, X_aux, Q, gamma, dir)[0];	  								 
Omega = 2*diagonal(-Q_matrix)';			  	
B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];	

obs_bloco = countc(y,X[][1]);	// conta o número de observações em cada intervalo [X_i, X_{i+1}). 
								// primeira linha obrigatoriamente é o número de observações em tempo inferior a w_i.
										
psi_aux = cond_Psi_MH(e, X ~ X_aux[][1],X[][2], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi, prioris_psi, par_priori_psi)[0];							

psi_aux2 = <>;
for(j = 0; j < sizeof(psi); ++j){	  
	psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];  // psi gerado pela condicional, será ele mesmo,
													   // caso não houve geração na respectiva condicional, o valor "anterior" será o último gerado
	psi_aux2 ~= psi_aux[j];	 //vetorizando psi, apenas para guardar os dados,
										 // mantendo os zeros na iteração onde não houve atualização
}




// inicio da calibração
					   
decl num_intervalos_vec, const_calib_vec, probs_aceit_vec, prob_limite, i2;

num_intervalos_vec = 0 | 10; //numero de intervalos para amostrar V
								// aqui, preciso colocar dois valores distintos
								// e que não podem ser possíveis para esse parâmetro em uma primeira rodada,
								// para poder entrar no "while"
const_calib_vec = zeros(1,columns(g_nao_constante)) | 2; // vetor das constantes de calibração do MH para psi
														 // (explicando: mesmo de cima!)

probs_aceit_vec = new array[0];	 // array das probs de aceitação para amostrar V
prob_limite = 0.25;				


// roda uma primeira vez: 

decl iterCalib_MCMC;

iterCalib_MCMC = calib_MCMC(calib, iterWVR, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2,valores_r, probs_r, classe_priori_r, hiperpar_r,
				prioris_psi, par_priori_psi);
i2 = 1;

X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][];


while(num_intervalos_vec[sizer(num_intervalos_vec)-1] .!=
			num_intervalos_vec[sizer(num_intervalos_vec)-2] ||
		sumr(const_calib_vec[sizer(const_calib_vec)-1][] .!=
			const_calib_vec[sizer(const_calib_vec)-2][]).!=0){	   

// calculo da matriz de covariancia empirica			
decl aux;
aux = cumulate(psi_tamanho' )';	
for(i1 = 0; i1 < columns(g_nao_constante); ++i1){
	matriz_cov[i1] = variance(iterCalib_MCMC[3][][aux[g_nao_constante[i1]] - psi_tamanho[g_nao_constante[i1]]:
								aux[g_nao_constante[i1]] - 1]) .* const_calib[i1];
	matriz_cov[i1] = mpd(matriz_cov[i1], 10^(-5));   // verificando as correlações empiricas, se são significativas (acima de 0.3)
}

// primeiro bloco: avaliar probs de aceitação de V e obter novo numero de intervalos.
aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));		

iterCalib_MCMC = calib_MCMC(calib, iterWVR, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2,valores_r, probs_r, classe_priori_r, hiperpar_r,
				prioris_psi, par_priori_psi);
i2 += 1;	   
				
X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][]; 

// avaliando probs de aceitação de V
// e mudando (ou não) o número de intervalos para amostrar V de acordo com a prob de aceitação
// (estratégia: utilizando a mediana das probs de aceitação em cada intervalo)

probs_aceit_vec |= {aceitCTMC ./ total_testesCTMC};		 

if(quantiler(aceitCTMC ./ total_testesCTMC) <= prob_limite){   
//	num_intervalos += 1;	   
}
			
// segundo bloco: obter probs de aceitação de Psi e obter nova constante de calibração.
aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));		

iterCalib_MCMC = calib_MCMC(calib, iterWVR, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, &aceitCTMC,
				&total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				&aceitPsi, &total_testesPsi, psi_aux2,valores_r, probs_r, classe_priori_r, hiperpar_r,
				prioris_psi, par_priori_psi);
i2 += 1;	
				
X = iterCalib_MCMC[0];	
Q_matrix = iterCalib_MCMC[1];
psi = iterCalib_MCMC[2];	
psi_aux2 = iterCalib_MCMC[3][calib-1][]; 
				 
const_calib = calibracao(const_calib, aceitPsi, total_testesPsi); // calcula a constante de calibração da matriz
																  // de cov de psi de acordo com a prob de aceitação de psi.
																  // Esta deve se manter entre 0.3 e 0.4.
																  
num_intervalos_vec |= num_intervalos;
const_calib_vec |= const_calib;	  
						  
}


// BURNIN

decl theta_post, psi_post, logposteriori, integral, logvero, tamanho_X, tempo_salto1;

theta_post = constant(.NaN, burnin+iter, e*e);
psi_post = constant(.NaN, burnin+iter, sumr(psi_tamanho));
logposteriori = constant(.NaN,burnin+iter,1);
integral = 	constant(.NaN,burnin+iter,1);
logvero = constant(.NaN, burnin+iter,1);
tamanho_X = constant(0,burnin+iter,1);
tempo_salto1 = constant(0,burnin+iter,1);

aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

for(i = 0; i < burnin; ++i){  	 

  for(j = 0; j < iterWVR; ++j){
	W = cond_W(X[][0:1], S, Q_matrix);
	intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];  // matriz em que cada linha são os intervalos entre saltos da ctmc(virtuais ou não)  [w_i, w_{i+1})						 
	obs_bloco = dropr(countc(y,W), rows(W));	  // conta o número de observações em cada intervalo [w_i, w_{i+1}). 
											  // primeira linha obrigatoriamente é zero.
																 				
	X = cond_VR(psi,num_intervalos,y,W, X, S,
		B_matrix, dist_inicial, &aceitCTMC, &total_testesCTMC,
		intervalos_tempos, obs_bloco,
		k, valores_r, probs_r, classe_priori_r, hiperpar_r);
  }

    if(rows(X) .== 1){ println("CTMC sem salto");
	tamanho_X[i] = 1;
	tempo_salto1[i] = S[1];
	}
	else{
  	tamanho_X[i] = rows(X);
	tempo_salto1[i] = X[1][1];
	}
	
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);
	Q_matrix = cond_theta(k, e, X_aux, Q_matrix, gamma, dir);
	logposteriori[i] = Q_matrix[1];  //adicionando a "parte 2 e 3" na logposteriori
	Q_matrix = Q_matrix[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];

	obs_bloco = countc(y,X[][1]);	// conta o número de observações em cada intervalo [X_i, X_{i+1}). 
									// primeira linha obrigatoriamente é o número de observações em tempo inferior a w_i.
				
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1],X[][2], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi,prioris_psi, par_priori_psi);

	logposteriori[i] +=	psi_aux[1]; //adicionando a "parte 4" na logposteriori
	psi_aux = psi_aux[0];

	psi_aux2 = <>;
	for(j = 0; j < sizeof(psi); ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];  // psi gerado pela condicional, será ele mesmo,
														   // caso não houve geração na respectiva condicional, o valor "anterior" será o último gerado
		psi_aux2 ~= psi_aux[j];  //vetorizando psi, apenas para guardar os dados,
										   // mantendo os zeros na iteração onde não houve atualização		
	}
	
	theta_post[i][] = vecr(Q_matrix)';
	psi_post[i][] = psi_aux2;
	integral[i] = sumc(logverossimilhanca(e,psi,y,obs_bloco,X,intervalos_tempos,range(0,rows(X)-1), X[][2])[1]);
	logposteriori[i] +=	logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, range (0, rows(X)-1), X[][2])[0]+
							log(dist_inicial[X_aux[0][0]]);	  // adicionando a "parte 1" da logpost
								 
	logvero[i] = logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, range (0, rows(X)-1),X[][2])[0];		
	decl aux5;						
	for(i1 = 0; i1 < rows(X); ++i1){	 // adicionando ultima parte da logposteriori, relativa a R
		aux5 = log(probs_r[X[i1][0]][][vecindex(valores_r[X[i1][0]] .== X[i1][2])]);
		if(aux5 == <>){
			logposteriori[i] += 0;
		}
		else{
			logposteriori[i] += aux5; // adicionando a "parte 5" da logpost (prob de R, caso discreto) 
		}
		if(classe_priori_r[X[i1][0]] == "R_gamma"){
			logposteriori[i] += (hiperpar_r[X[i1][0]][0][0] - 1)
											*log(X[i1][2]) - hiperpar_r[X[i1][0]][0][1]*X[i1][2];
						   			
		}
	}
	
}	
	
// INICIO DA CADEIA!


decl X_iter, lambda_iter, int_cred;

lambda_iter = constant(0, 1, columns(grid));  
int_cred = constant(0,iter,stepgrid+1);
	
aceitCTMC = zeros(1,num_intervalos);
total_testesCTMC = zeros(1,num_intervalos);

aceitPsi = zeros(1, columns(g_nao_constante));	  
total_testesPsi = zeros(1, columns(g_nao_constante));


decl lambda_previsao, PP_previsao, integral_previsao;
lambda_previsao = constant(0, 1, columns(range(S[1],S[1]+S_pred,S_pred./100)));
PP_previsao = constant(0,iter,1);
integral_previsao = constant(0,iter,1);

for(i = 0; i < iter; ++i){	  	

  for(j = 0; j < iterWVR; ++j){
	W = cond_W(X[][0:1], S, Q_matrix);
	intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];  // matriz em que cada linha são os intervalos entre saltos da ctmc(virtuais ou não)  [w_i, w_{i+1})						 
	obs_bloco = dropr(countc(y,W), rows(W));	  // conta o número de observações em cada intervalo [w_i, w_{i+1}). 
											  // primeira linha obrigatoriamente é zero.
																 				
	X = cond_VR(psi,num_intervalos,y,W, X,S,
			 B_matrix,dist_inicial,&aceitCTMC,&total_testesCTMC,
			intervalos_tempos,obs_bloco,
			k,valores_r,probs_r,classe_priori_r,hiperpar_r);	 
  }

  	if(rows(X) .== 1){ 
	tamanho_X[burnin+i] = 1;
	tempo_salto1[burnin+i] = S[1];
	}
	else{
	tamanho_X[burnin+i] = rows(X);
	tempo_salto1[burnin+i] = X[1][1];
	}
			
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);
	Q_matrix = cond_theta(k, e, X_aux, Q_matrix, gamma, dir);
	logposteriori[burnin+i] = Q_matrix[1];	  //adicionando a "parte 2 e 3" na logposteriori
	Q_matrix = Q_matrix[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];

	obs_bloco = countc(y,X[][1]);	// conta o número de observações em cada intervalo [X_i, X_{i+1}). 
									// primeira linha obrigatoriamente é o número de observações em tempo inferior a w_i.
				
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1],X[][2], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi, &total_testesPsi,prioris_psi, par_priori_psi);
	logposteriori[burnin+i] += psi_aux[1];
	psi_aux = psi_aux[0];
	
	psi_aux2 = <>;
	for(j = 0; j < sizeof(psi); ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];  // psi gerado pela condicional, será ele mesmo,
														   // caso não houve geração na respectiva condicional, o valor "anterior" será o último gerado
		psi_aux2 ~= psi_aux[j];  //vetorizando psi, apenas para guardar os dados,
										   // mantendo os zeros na iteração onde não houve atualização		
	}
	
	// guardando resultados!!!
	
	decl aux1, aux2, aux3, aux4, aux5, lambda;
	aux1 = sumc(X[][1] .<= grid)-1;										  
	aux3 = <>;
	aux4 = <>;
	for(j = 0;j<e;++j){
		aux3 ~= X[aux1][1];
		aux4 ~= X[aux1][2];
	}
			
	lambda = func_g(psi, grid',aux3,aux4);	 
	X_iter = X[aux1][0]';	
	aux2 = selectrc(lambda,range(0,columns(grid)-1, 1)', X_iter);  
	lambda_iter += aux2;
	int_cred[i][] = aux2[range(0,stepgrid)];
							
	theta_post[burnin+i][] = vecr(Q_matrix)';
	psi_post[burnin+i][] = psi_aux2;

	integral[burnin+i] = sumc(logverossimilhanca(e,psi,y,obs_bloco,X,intervalos_tempos,range(0,rows(X)-1), X[][2])[1]);
	logposteriori[burnin+i] +=	logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, range (0, rows(X)-1), X[][2])[0]+
							log(dist_inicial[X_aux[0][0]]);	 
	logvero[burnin+i] = logverossimilhanca(e,psi,y,obs_bloco,X,intervalos_tempos,range(0,rows(X)-1),X[][2])[0];						   																														
							
	for(i1 = 0; i1 < rows(X); ++i1){	 // adicionando ultima parte da logposteriori, relativa a R
		aux5 = log(probs_r[X[i1][0]][][vecindex(valores_r[X[i1][0]] .== X[i1][2])]);
		if(aux5 == <>){
			logposteriori[burnin+i] += 0;
		}
		else{
			logposteriori[burnin+i] += aux5;  // adicionando a "parte 5" da logpost (prob de R, caso discreto) 
		}
		if(classe_priori_r[X[i1][0]] == "R_gamma"){
			logposteriori[burnin+i] += (hiperpar_r[X[i1][0]][0][0] - 1)
											*log(X[i1][2]) - hiperpar_r[X[i1][0]][0][1]*X[i1][2];
		}
	}	   
	
	if(logposteriori[burnin+i] == .NaN){println(logposteriori[burnin+i]);
				println(avaliar_lambda(X | (X[rows(X)-1][0] ~ S[1] ~ X[rows(X)-1][2]), psi, e, X[][2]));}
													
	if(fmod(i+1,10000)==0){println("time for ", burnin+i+1," iter: ",timespan(time1));}

	if(S_pred > 0){
	decl previsao;

	previsao = gerar_processo(S, S_pred, E, Q_matrix, psi, e, valores_r, probs_r,
				 	range(S[1],S[1]+S_pred,S_pred./100),X[rows(X)-1][]);
	PP_previsao[i][] = previsao[1];
	integral_previsao[i][] = previsao[2];
	}	
	
}							


	println("accept probability (V_0, V):", aceitCTMC ./ total_testesCTMC);
	println("accept probability Psi:", aceitPsi ./ total_testesPsi);

//	savemat("lambda_mean.mat", lambda_iter ./ iter, 1);
//	savemat("lambda_ci.mat", quantilec(int_cred,<0.025,0.975>), 1);
//	savemat("psi.mat",psi_post,1);
//	savemat("logposterior.mat",logposteriori,1);
	savemat("loglike.mat",logvero,1);
//	savemat("integral.mat",integral,1);
//	savemat("CTMCjumps.mat",tamanho_X,1);
//	savemat("firstjump.mat",tempo_salto1,1);
//	
//	if(rows(gamma) .>= 1){
//	savemat("theta.mat",theta_post,1);
//	}
//	
//	if(S_pred > 0){	 
//	savemat("Npredic.mat", PP_previsao, 1);	
//	savemat("integralpredic.mat", integral_previsao, 1);
//	}

}
