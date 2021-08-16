#include <functional_forms.ox>


mpd(const m, const tol){
  decl esv,esve,delta,tau,dm;

  eigensym(m,&esv,&esve);

  delta=2*tol;
  
  tau=((delta-esv) .<0 .?0 .:(delta-esv));

  dm=esve*diag(tau)*esve';
  

  return (m+dm);
}


rmvnorm(const r,const c ,const mu, const sigma){
  decl u, A, B;

  B = (sigma+sigma')/2;		 
  A = choleski(B);	 
  u = mu + rann(r,c)*A;

  return u;		
 }

  


cond_W(const X, const S, const Q){

decl e;
e = columns(Q);
 								
decl tempo_espera, Omega, intensity, u_aux, U, i;

tempo_espera = (X | S[1])[1:][1] - X[][1];
Omega = 5*diagonal(-Q)'; 			 
intensity = Omega[X[][0]] + minr(Q[X[][0]][]); 	   
						
u_aux = ranpoisson(rows(intensity), 1, intensity .* tempo_espera);	
U = <>;									
for(i = 0; i < rows(intensity); ++i){
	U = U | (X[i][1] + ranu(u_aux[i],1) .* tempo_espera[i]); 
}			

decl W;
W = sortc(U | X[][1]);	
						
return W;		  

}

   

cond_V(const psi,const y, const W, const B_matrix, const S, const b2){

decl V_aux, probs_V, traj, prob_trans, indices, obs_bloco, integral, soma_loglambda, j;
V_aux = setupper(zeros(rows(W)-1,rows(W)),1);	 
probs_V = constant(.NaN, 1, rows(W)-1); 

for(j = 0; j < rows(W)-1; ++j){
   	traj = V_aux[j][]' ~ W;	 
	decl jj;
	prob_trans = 0; 
	for(jj = 0; jj < rows(W)-1; ++jj){
		prob_trans += log(B_matrix[traj[jj][0]][traj[jj+1][0]]); 
	}	
	
	indices = 1 | (traj[1:(rows(traj)-1)][0] .== traj[0:(rows(traj)-2)][0] .? 0 .: 1);	
	traj = traj[vecindex(indices .== 1)][];	
	
	obs_bloco = countc(y,traj[][1]); 
	integral = sumr(int_g(psi,traj[1][1],S, b2));

	if(obs_bloco[1] .== 0){
		soma_loglambda = 0;
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[0]:cumulate(obs_bloco)[1]-1], traj[1][1],b2)))[0];
	}
	if(obs_bloco[2] .== 0){
		soma_loglambda += 0;
	}
	else{
		soma_loglambda += sumc(log(func_g(psi,y[cumulate(obs_bloco)[1]:cumulate(obs_bloco)[2]-1], traj[1][1],b2)))[1];
	}
	
	probs_V[j] = -integral + soma_loglambda + prob_trans;

	if(func_g(psi,traj[1][1],traj[1][1],b2)[0] .< b2){	
		probs_V[j] = 1;
	}
	  	if(probs_V == .NaN){println(probs_V);}
}
 
decl aux;
aux = maxr(probs_V);
probs_V = exp(probs_V - log(sumr(exp(probs_V - aux))) - aux);
					if(probs_V == .NaN){println(probs_V);}
decl V;
V = V_aux[vecindex(ranmultinomial(1,probs_V) .== 1)][]' ~ W; 
indices = 1 | (V[1:(rows(V)-1)][0] .== V[0:(rows(V)-2)][0] .? 0 .: 1);	
V = V[vecindex(indices .== 1)][];

return V;

}	


cond_theta(const k, const e, const X_aux, const Q){

decl theta1,j;
theta1 = diagonal(-Q);

decl aux3, mj, deltaj, Aj, Bj, aux4, aux1log, aux2log;

aux1log = constant(.NaN, 1, k); 
aux2log = constant(.NaN, 1, k); 

for(j = 0; j < k; ++j){					 
				 
aux3 = X_aux[vecindex(X_aux[][0] .== j .|| X_aux[][0] .== j + k)][];  																	 
mj = rows(aux3); 
deltaj = sumc(aux3)[1];  

if(X_aux[rows(X_aux)-1][0] .== j .|| X_aux[rows(X_aux)-1][0] .== j + k){ 
aux1log[j] = mj - 1;  
}
else{ 	  
aux1log[j] = mj;  
}

aux2log[j] = deltaj; 

}

decl logP2;	   

logP2 = sumr(log(theta1[0:(k-1)]) .* aux1log) - sumr(theta1[0:(k-1)] .* aux2log);

decl i1, i2, aux1, mdots;
if(rows(X_aux) > 1){   
	mdots = zeros(2*k, 2*k);
	aux1 = X_aux[:(rows(X_aux)-2)][0] ~ X_aux[1:][0];	
	for(i1 = 0; i1 < e; ++i1){	
		for(i2 = 0; i2 < e; ++i2){
			decl x;
			x = rows(vecindex((aux1[][0] .== i1) .&& (aux1[][1] .== i2)));
			mdots[i1][i2] = x;
		}
	}	 
}
else{	
	mdots = zeros(2*k, 2*k);
}

decl mdots_aux;
mdots_aux = mdots[][0:(k-1)] + mdots[][k:(2*k-1)]; 
mdots_aux = mdots_aux[0:(k-1)][] + mdots_aux[k:(2*k-1)][]; 
				
decl logP3, aux3log;
aux3log = setdiagonal(Q ./ (diagonal(-Q)'),0);			 

if(rows(aux3log) != 2*k){ 	 
aux3log ~= constant(0,rows(aux3log),2*k-rows(aux3log));
}

aux3log = aux3log[][0:k-1] + aux3log[][k:2*k-1];
aux3log = aux3log[0:k-1][]; 		

aux3log = mdots_aux  .* log(aux3log);	 

logP3 = 0;
for(j = 0; j < k; ++j){
logP3 += sumr(deletec(aux3log[j][]));
}				 
								  
return logP2+logP3;

}



logverossimilhanca(psi,y, obs_bloco, X, S, b2){	  

decl integral, soma_loglambda, logvero;

integral = sumr(int_g(psi,X[1][1],S, b2));

if(obs_bloco[1] .== 0){
	soma_loglambda = 0;
}
else{
	soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[0]:cumulate(obs_bloco)[1]-1], X[1][1],b2)))[0];
}
if(obs_bloco[2] .== 0){
	soma_loglambda += 0;
}
else{
	soma_loglambda += sumc(log(func_g(psi,y[cumulate(obs_bloco)[1]:cumulate(obs_bloco)[2]-1], X[1][1],b2)))[1];
}
	
logvero = -integral + soma_loglambda; 
			  
return {logvero, integral};	  

}


teste_mhPsi(const psi, const logvero, const psi_new, const logvero_new, const aceit, const total_testes, const j){

total_testes[0][j] += 1;		

decl unif, alpha;
unif = ranu(1,1);
alpha = min(0,logvero_new - logvero);
			
if(log(unif) .< alpha){	 
	aceit[0][j] += 1;	 
	return psi_new;
}
else{	  
	return psi; 
}

}




dens_priori_psi(const psi, const prioris_psi, const par_priori_psi){

decl psi_vec, i;

psi_vec = psi[0];
for(i = 1; i < rows(psi); ++i){
	psi_vec ~= psi[i];
}				
								 
decl logdens_priori;
logdens_priori = 0;
for(i = 0; i < rows(prioris_psi);++i){
	if(prioris_psi[i] == "normal"){
		logdens_priori += -1/(2*par_priori_psi[i][0][1])*(psi_vec[0][i]-par_priori_psi[i][0][0])^2;
	}
	if(prioris_psi[i] == "uniform"){
		logdens_priori += 0;
	}
}
return logdens_priori;
}



cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, aceit, total_testes,
					prioris_psi, par_priori_psi, b2){

decl j, ind, psi_new;	   
decl dens_transicao, dens_transicao_new;		  

for(j = 0; j < k; ++j){
			 
		psi_new = psi;								 
		psi_new[j] = rmvnorm(1, columns(psi[j]), psi[j], matriz_cov[j]);  
									
		if(condicao_sinal_psi(psi_new) .== 1){																
			dens_transicao = logverossimilhanca(psi,y, obs_bloco, X, S, b2)[0] 	
							 + dens_priori_psi(psi, prioris_psi, par_priori_psi);
			dens_transicao_new = logverossimilhanca(psi_new,y, obs_bloco, X, S, b2)[0]	
					   		 + dens_priori_psi(psi_new, prioris_psi, par_priori_psi);							 
			psi[j] = teste_mhPsi(psi[j], dens_transicao, psi_new[j], dens_transicao_new, &aceit[0], &total_testes[0], j);
		
		}
		else{ 	   			  
			total_testes[0][j] += 1;	
		}

}

decl logP4;
logP4 = dens_priori_psi(psi, prioris_psi, par_priori_psi);

return {psi,logP4};

}



calib_MCMC(calib, X, S, Q, psi, y, Omega, B_matrix, k, e, psi_tamanho, matriz_cov,
				aceitPsi, total_testesPsi, prioris_psi, par_priori_psi, trunc_W, b2){

decl psi_calib;				
psi_calib = constant(.NaN, calib, sumr(psi_tamanho));
															
decl i, W, X_aux, obs_bloco, psi_aux, j, psi_aux2;

for(i = 0; i < calib; ++i){	  

W = cond_W(X, S, Q);
W = W[0][] | (W[vecindex(W .>=trunc_W[0] .&& W .<=trunc_W[1])][]);
X = cond_V(psi, y, W, B_matrix, S, b2);
X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]); 
obs_bloco = countc(y,X[][1]);
psi_aux = cond_Psi_MH(X, S, k, psi, matriz_cov, y, obs_bloco, &aceitPsi[0], &total_testesPsi[0],
						prioris_psi, par_priori_psi, b2)[0];							
psi_aux2 = <>;
for(j = 0; j < k; ++j){
	psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];  
	psi_aux2 ~= psi_aux[j];	 
}
	  
	psi_calib[i][] = psi_aux2;
	
}	  

return {X, psi, psi_calib};

}


   
calibracao(const_calib, const aceitPsi, const total_testesPsi){
						
decl k, prob_aceit;
					   
	prob_aceit = aceitPsi[0] ./ total_testesPsi[0];	 
	const_calib[0][0] = prob_aceit .< 0.15 .? const_calib[0][0] .* 0.9 .:
						prob_aceit .> 0.3 .? const_calib[0][0] .* 1.1 .: const_calib[0][0];
								
	prob_aceit = aceitPsi[1] ./ total_testesPsi[1];	 
	const_calib[0][1] = prob_aceit .< 0.2 .? const_calib[0][1] .* 0.9 .:
						prob_aceit .> 0.35 .? const_calib[0][1] .* 1.1 .: const_calib[0][1];
								
return const_calib;

}



predicao(psi,T,S, lambda_prev,b2){

decl gamma_g, func_inversa;
							 
gamma_g	= (psi[0][1] + psi[0][0]*probn(psi[0][3] + sqrt(2*M_PI)/psi[0][0] * psi[0][2] * T) - b2)./(probn(psi[1][1]));
func_inversa = -(quann((lambda_prev-b2)/gamma_g)-psi[1][1]) ./ (psi[1][0].*sqrt(2*M_PI)./gamma_g) + T;

decl tam_grid_int, grid_int, int2;
tam_grid_int = 100;
grid_int = range(S[1],func_inversa,(func_inversa-S[1])/tam_grid_int);	
int2 = 	func_g(psi,meanr(grid_int[0:(tam_grid_int - 1)]' ~ grid_int[1:tam_grid_int]'),T,b2)[][1];
int2 = sumc(int2 .* ((func_inversa-S[1])/tam_grid_int));

return {func_inversa, int2};

}