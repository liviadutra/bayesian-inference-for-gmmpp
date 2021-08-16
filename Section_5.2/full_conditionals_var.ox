#include <functional_forms.ox>


mpd(const m, const tol){
  decl esv,esve,delta,tau,dm;

  eigensym(m,&esv,&esve);

  delta=2*tol;
  
  tau=((delta-esv) .<0 .?0 .:(delta-esv));

  dm=esve*diag(tau)*esve';
  

  return (m+dm);
}



avaliar_lambda(X, psi, e, R){	  

	decl aux, epsilon, delta, i, inicial, final;
	aux = rows(X)-2;
	epsilon = 10^(-10); 

	delta = X[][1];
	for(i = 1; i < e; ++i){
		delta ~= X[][1];
	}
	delta = delta[:aux][]; 
	
	inicial = func_g(psi, X[:aux][1], delta, R);	 
	inicial = selectrc(inicial,range(0,aux),X[:aux][0])'; 
	
	final = func_g(psi, X[1:][1]-epsilon, delta, R);	
	final = selectrc(final,range(0,aux),X[:aux][0])'; 	 

	return inicial ~ final;		 

}

  
rmvnorm(const r,const c ,const mu, const sigma){
  decl u, A, B;

  B = (sigma+sigma')/2;		 
  
  A = choleski(B);	 
		   
  u = mu + rann(r,c)*A;

  return u;		
}


ranexpt(const lambda, const b){
  decl x;

  x = -(1 ./ lambda) .* log(1 - ranu(1,1).*(1 - exp(-lambda .* b)) );
  
  return x;			
} 


R_discreta(r, vec_valores, vec_prob){
decl aux;
aux = vecindex(vec_valores .== r);
return vec_prob[aux];
}

 
logConst_R_discreta(tam_r,  psi, delta, valores_r, probs_r, j, intervalos, obs_bloco,y,hiperpar_r){

decl aux3, num_r, avaliar_func, integral, soma_loglambda;
aux3 = constant(.NaN,1,tam_r);	 

for(num_r = 0; num_r < tam_r; ++num_r){	
avaliar_func = (func_g(psi, intervalos[0], delta, valores_r[num_r])[j] .< 0) .||
						(func_g(psi, intervalos[1], delta, valores_r[num_r])[j] .< 0);	
			
integral = int_g(psi,intervalos, delta, valores_r[num_r])[j];  
				
if(obs_bloco .== 0){
	soma_loglambda = <0>;	
}
else{	
	soma_loglambda = sumc(log(func_g(psi,y,delta, valores_r[num_r])))[j];
}	  
	
if(avaliar_func == 1){ 
	aux3[num_r] = -.Inf;
	}
else{								
	aux3[num_r] = -integral + soma_loglambda + log(probs_r[num_r]);
	}

}		  

return aux3;

}

  

logConst_R_uniforme(tam_r,  psi, delta, valores_r, probs_r, j, intervalos, obs_bloco,y,hiperpar_r){
	decl alfa, beta, log_const_gamma, i, log_sum;
	alfa = obs_bloco + 1;
	beta = intervalos[1] - intervalos[0];
	
	log_sum = 0;
	for(i = 1; i < alfa; ++i){
	log_sum += log(i);
	}	
	log_const_gamma = matrix (alfa * log(beta) - log_sum);

	return log_const_gamma;
}

logConst_R_gamma(tam_r,  psi, delta, valores_r, probs_r, j, intervalos, obs_bloco,y,hiperpar_r){
	decl alfa, beta, log_const_gamma, i, log_sum;
	alfa = obs_bloco + hiperpar_r[0];
	beta = intervalos[1] - intervalos[0] + hiperpar_r[1];
	
	log_sum = 0;
	for(i = 1; i < alfa; ++i){
	log_sum += log(i);
	}	
	log_const_gamma = matrix (alfa * log(beta) - log_sum);

	return log_const_gamma;
}

		

gerando_R_discreta(constantes_auxV, V, valores_r, intervalos, obs_bloco, hiperpar_r){

	decl prob_R, r;
	prob_R = exp(constantes_auxV[V] - log(sumr(exp(constantes_auxV[V] - maxr(constantes_auxV[V]))))
								- maxr(constantes_auxV[V]));

	r = valores_r[vecindex(ranmultinomial(1,prob_R) .== 1)];

	return r;
}


gerando_R_uniforme(constantes_auxV, V, valores_r, intervalos, obs_bloco, hiperpar_r){

	decl r;										   

	r = rangamma(1,1,obs_bloco+1, intervalos[1] - intervalos[0]);

	return r;
}


gerando_R_gamma(constantes_auxV, V, valores_r, intervalos, obs_bloco, hiperpar_r){

	decl r;										   

	r = rangamma(1,1,obs_bloco + hiperpar_r[0],	intervalos[1] - intervalos[0] + hiperpar_r[1]);;

	return r;
}

  

gerando_vr(const intervalos_tempos, const obs_bloco,
			const psi,const y, const W, const B_matrix, const dist_inicial, const e, const k,
			const valores_r, const probs_r, const classe_priori_r, const hiperpar_r){

decl V, integral, avaliar_func,
		soma_loglambda, i, delta, j, constantes_new, R_aux;
								
V = constant(.NaN, rows(W)-1, 1);   // V trajetoria da DTMC subordinada aos tempos W.
R_aux = constant(.NaN, rows(W)-1, 1);

i = 0;	
delta = constant(W[i], 1, e);

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_Final;							
constantes_Final = 0;		

decl constantes_auxV, constantes_V, aux4, constantes_auxFinal, prob_V;
constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

for(j = 0; j < k; ++j){		
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 				 
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(dist_inicial[j]) + constantes_V[j]; 
}

constantes_Final += 0; 
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);
	 
R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
									intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
									
decl delta_aux;
delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

for(i = 1; i < rows(W)-1; ++i){

constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
	
for(j = 0; j < k; ++j){

	if(V[i-1] .== j){		
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){						
				constantes_auxFinal[j] = -.Inf;
			}														
}

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 																		   
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								
V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

}

V = V ~ W[0:rows(W)-2] ~ R_aux;	
V = deleter(V);

return {V, constantes_Final};

}


 
gerando_vr_prim(const intervalos_tempos, const obs_bloco, const V_b_cond,
			const psi,const y, const W, const B_matrix, const dist_inicial, const e, const k,
			const valores_r, const probs_r, const classe_priori_r, const hiperpar_r){

decl V, integral, avaliar_func,
		soma_loglambda, i, delta, j, constantes_new, R_aux;
										
V = constant(.NaN, rows(W), 1);  
R_aux = constant(.NaN, rows(W), 1);
		
V[rows(W)-1][0] = V_b_cond; 
	  
i = 0;	 
delta = constant(W[i], 1, e);

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_Final;							
constantes_Final = 0;		

decl constantes_auxV, constantes_V, aux4, constantes_auxFinal, prob_V;
constantes_auxV = new array[k]; 	 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

for(j = 0; j < k; ++j){	 			
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 				
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(dist_inicial[j]) + constantes_V[j];
}

constantes_Final +=	0; 
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);
	 
R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
									intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
									
decl delta_aux;
delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;		 
			
for(i = 1; i < rows(W)-2; ++i){

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

for(j = 0; j < k; ++j){	 
	if(V[i-1] .== j){	
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
			   		
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 	
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){					
				constantes_auxFinal[j] = -.Inf;
			}														
																	
}

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 
																		   
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;		 

}

i = rows(W)-2;

if(i .== 0){ }

else{

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

for(j = 0; j < k; ++j){	 
	if(V[i-1] .== j){		   
	decl aux1, aux2;
	aux1 = unit(1,1);	
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{	  
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 	 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + log(B_matrix[j][V[i+1]]')
										+ constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
													
}

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 
																		   
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								
prob_V[V[rows(W)-1]] = 0;	
V[i] = 	vecindex(ranmultinomial(1,prob_V ./ sumr(prob_V)) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

}


V = V ~ W ~ R_aux;	 
V = deleter(V);

return {V, constantes_Final};

}


 
gerando_vr_inter(const intervalos_tempos, const obs_bloco, const V_b_cond,	delta,
			const psi,const y, const W, const B_matrix, const dist_inicial, const e, const k,
			const valores_r, const probs_r, const classe_priori_r, const hiperpar_r){		 
					 							  
decl V, integral, avaliar_func,
		soma_loglambda, i, j, constantes_new, R_aux;
										
V = constant(.NaN, rows(W), 1);
R_aux = constant(.NaN, rows(W), 1);
		
V[rows(W)-1][0] = V_b_cond[1];

i = 0;
decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_Final;							
constantes_Final = 0;		

decl constantes_auxV, constantes_V, aux4, constantes_auxFinal, prob_V;
constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);
	 
for(j = 0; j < k; ++j){

	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);			
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V_b_cond[0]][j]) + constantes_V[j]; 
}

constantes_Final +=	0; 
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);
	 
R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
									intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
								
decl delta_aux;
delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;		

for(i = 1; i < rows(W)-2; ++i){

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
	
for(j = 0; j < k; ++j){	
	if(V[i-1] .== j){		
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
			   		
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j];
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
																
}

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 
																		   
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;		 

}
 
i = rows(W)-2;

if(i .== 0){ }

else{

constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
	
for(j = 0; j < k; ++j){	 
	if(V[i-1] .== j){		 
	decl aux1, aux2;
	aux1 = unit(1,1);	
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{	  
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 	 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + log(B_matrix[j][V[i+1]]')
										+ constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){		
				constantes_auxFinal[j] = -.Inf;
			}														
												
}

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 			   
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								
prob_V[V[rows(W)-1]] = 0;

V[i] = 	vecindex(ranmultinomial(1,prob_V ./ sumr(prob_V)) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

}


V = V ~ W ~ R_aux;	 
V = deleter(V);

return {V, constantes_Final};

}



gerando_vr_ultimo(const intervalos_tempos, const obs_bloco, const V_b_cond,	const R_b_cond, delta,
			const psi,const y, W, const B_matrix, const dist_inicial, const e, const k,
			const valores_r, const probs_r, const classe_priori_r, const hiperpar_r){		 
			
decl W_aux, V, integral, avaliar_func,
		soma_loglambda, i, j, constantes_new, R_aux;

W_aux = W[0]; 	
V = constant(.NaN, rows(W)-1, 1);   
R_aux = constant(.NaN, rows(W)-1, 1); 
	
i = 0;

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_Final;							
constantes_Final = 0;		

decl constantes_auxV, constantes_V, aux4, constantes_auxFinal, prob_V;
constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

	
for(j = 0; j < k; ++j){

	if(V_b_cond .== j){ 		
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = R_b_cond; 	
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V_b_cond][j]) + constantes_V[j]; 
}	

constantes_Final +=	0; 
					
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);
	 
if(V[i] .!= V_b_cond){	
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}	
else{
	R_aux[i] = R_b_cond;
	W[i] = delta[vecindex(delta .!= intervalos_tempos[i][0])];
}

decl delta_aux;
delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;		

for(i = 1; i < rows(W)-1; ++i){

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{	
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

for(j = 0; j < k; ++j){	
	if(V[i-1] .== j){ 	
	decl aux1, aux2;
	aux1 = unit(1,1);

	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{					
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}	  
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
																	
}	  

constantes_Final += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal);
prob_V = exp(constantes_auxFinal - log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) -
								maxr(constantes_auxFinal));								

V[i] = 	vecindex(ranmultinomial(1,prob_V) .== 1);

if(V[i] .!= V[i-1]){
	R_aux[i] = sprint("gerando_",classe_priori_r[V[i]])(constantes_auxV, V[i],valores_r[V[i]],
										intervalos_tempos[i][], obs_bloco[i+1], hiperpar_r[V[i]]);
}

delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

}				

W[0] = W_aux;
V = V ~ dropr(W, (rows(W)-1)) ~ R_aux;	
V = deleter(V);	 
							
return {V, constantes_Final};

}




constantes_func_cond(const intervalos_tempos, const obs_bloco, delta,
				const psi,const y, const W, const X, const R, const B_matrix, const dist_inicial,
				const e, const k, const valores_r, const probs_r, const classe_priori_r,
				const hiperpar_r){	   

decl V, R_aux, constantes, i, j;

V = constant(.NaN, rows(W)-1, 1); 		 				 
R_aux = constant(.NaN, rows(W)-1, 1);
constantes = 0;

i = 0;

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_auxV, constantes_V, aux4, aux5, constantes_auxFinal;
constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

	
for(j = 0; j < k; ++j){	 			
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]); 
		 				
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(dist_inicial[j]) + constantes_V[j]; 
}

constantes += 0; 
aux5 = maxc(vecindex(W[i] .>= X[][1]));				
V[i] = X[aux5][0]; 
R_aux[i] = R[aux5];

decl delta_aux;
delta_aux = delta[V[i]];		 		  			
delta = constant(W[i+1], 1, e);					
delta[V[i]] = delta_aux;			  				

for(i = 1; i < rows(W)-2; ++i){

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
 
for(j = 0; j < k; ++j){

	if(V[i-1] .== j){	
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){					
				constantes_auxFinal[j] = -.Inf;
			}														
																	
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 

aux5 = maxc(vecindex(W[i] .>= X[][1]));								
V[i] = X[aux5][0];
R_aux[i] = R[aux5];


delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

}

i = rows(W)-2;

if(i .== 0){ }

else{

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

for(j = 0; j < k; ++j){	  

	if(V[i-1] .== j){	
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 		
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) +	log(B_matrix[j][X[rows(X)-1][0]]') +
											constantes_V[j];
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
														
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 

}

return constantes;

}





constantes_func(const intervalos_tempos, const obs_bloco, delta,
				const psi,const y, const W, const X, const R, const B_matrix, const dist_inicial,
				const e, const k, const valores_r, const probs_r, const classe_priori_r,
				const hiperpar_r){	   

decl V, R_aux, constantes, i, j;

V = constant(.NaN, rows(W)-1, 1);  		 				 
R_aux = constant(.NaN, rows(W)-1, 1);
constantes = 0;

i = 0;

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_auxV, constantes_V, aux4, aux5, constantes_auxFinal;
constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);


for(j = 0; j < k; ++j){	 
						 			
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]); 
		 				
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(dist_inicial[j]) + constantes_V[j];
}

constantes += 0;
aux5 = maxc(vecindex(W[i] .>= X[][1]));				
V[i] = X[aux5][0]; 
R_aux[i] = R[aux5];

decl delta_aux;
delta_aux = delta[V[i]];		 		  			
delta = constant(W[i+1], 1, e);					
delta[V[i]] = delta_aux;			  				

for(i = 1; i < rows(W)-2; ++i){

constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
   
for(j = 0; j < k; ++j){

	if(V[i-1] .== j){		 
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j];
			if(constantes_auxFinal[j] .== .NaN){					
				constantes_auxFinal[j] = -.Inf;
			}														
																	
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal);

aux5 = maxc(vecindex(W[i] .>= X[][1]));								
V[i] = X[aux5][0];
R_aux[i] = R[aux5];


delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

}

i = rows(W)-2;

if(i .== 0){ }

else{

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
			

for(j = 0; j < k; ++j){	 

	if(V[i-1] .== j){		
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
																
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 

}

return constantes;

}






constantes_func_cond_ultimo(const intervalos_tempos, const obs_bloco, V_b_cond, R_b_cond, delta,
				const psi,const y, W, const X, const R, const B_matrix, const dist_inicial,
				const e, const k, const valores_r, const probs_r, const classe_priori_r,
				const hiperpar_r){	   

decl W_aux, V, R_aux, constantes, i, j;

W_aux = W[0];
V = constant(.NaN, rows(W)-1, 1);  			 				 
R_aux = constant(.NaN, rows(W)-1, 1);
constantes = 0;

i = 0;

decl y_aux;
if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

decl constantes_auxV, constantes_V, aux4, aux5, constantes_auxFinal;
constantes_auxV = new array[k]; 
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

for(j = 0; j < k; ++j){

	if(V_b_cond .== j){		
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = R_b_cond;
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]); 
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(dist_inicial[j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){					
				constantes_auxFinal[j] = -.Inf;
			}														
																	
}

constantes += 0; 
aux5 = maxc(vecindex(W[i] .>= X[][1]));			
V[i] = X[aux5][0];

if(V[i] .!= V_b_cond){
	R_aux[i] = R[aux5];	 
}	
else{					
	V[i] = V_b_cond;
	R_aux[i] = R_b_cond;	
	W[i] = delta[vecindex(delta .!= intervalos_tempos[i][0])];
}


decl delta_aux;
delta_aux = delta[V[i]];		 		  			
delta = constant(W[i+1], 1, e);					
delta[V[i]] = delta_aux;			  				

for(i = 1; i < rows(W)-2; ++i){

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}
		 
for(j = 0; j < k; ++j){	 
	if(V[i-1] .== j){	
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){				
				constantes_auxFinal[j] = -.Inf;
			}														
																
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 

aux5 = maxc(vecindex(W[i] .>= X[][1]));								
V[i] = X[aux5][0];
R_aux[i] = R[aux5];


delta_aux = delta[V[i]];	
delta = constant(W[i+1], 1, e);	  
delta[V[i]] = delta_aux;	

}

i = rows(W)-2;

if(i .== 0){ }

else{

constantes_auxV = new array[k];
constantes_V = constant(.NaN,1,k);
constantes_auxFinal = constant(.NaN,1,k);

if(obs_bloco[i+1] .== 0){
	y_aux = <>;
}
else{
	y_aux = y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1];
}

for(j = 0; j < k; ++j){

	if(V[i-1] .== j){	
	decl aux1, aux2;
	aux1 = unit(1,1);
	aux2 = matrix(deleter(R_aux)[rows(deleter(R_aux))-1]);
	constantes_auxV[j] = logConst_R_discreta
							(1,  psi, delta, aux2,
							aux1, j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
	}
	else{
	constantes_auxV[j] = sprint("logConst_",classe_priori_r[j])
							(rows(valores_r[j]),  psi, delta, valores_r[j],
							probs_r[j], j, intervalos_tempos[i][], obs_bloco[i+1],
							y_aux,hiperpar_r[j]);
		 			
	}
	aux4 = maxr(constantes_auxV[j]); 

	constantes_V[0][j] = log(sumr(exp(constantes_auxV[j] - aux4))) + aux4;			

	constantes_auxFinal[j] = log(B_matrix[V[i-1]][j]) + constantes_V[j]; 
			if(constantes_auxFinal[j] .== .NaN){					
				constantes_auxFinal[j] = -.Inf;
			}														
																
}

constantes += log(sumr(exp(constantes_auxFinal - maxr(constantes_auxFinal)))) +
								maxr(constantes_auxFinal); 

}

return constantes;

}





cond_W(const X, const S, const Q){

decl e;
e = columns(Q);
									
decl tempo_espera, Omega, intensity, u_aux, U, i;

tempo_espera = (X | S[1])[1:][1] - X[][1];
Omega = 2*diagonal(-Q)'; 			 
intensity = Omega[X[][0]] + minr(Q[X[][0]][]); 
						
u_aux = ranpoisson(rows(intensity), 1, intensity .* tempo_espera);
U = <>;									
for(i = 0; i < rows(intensity); ++i){
	U = U | (X[i][1] + ranu(u_aux[i],1) .* tempo_espera[i]);
}			

decl W;
W = sortc(U | X[][1]) | S[1];

return W;		  

}


				 



teste_mh(const X, const constantes, const X_new, const constantes_new, aceit, total_testes, const i){
	
total_testes[0][i] += 1;

decl unif, alpha;
unif = ranu(1,1);
alpha = min(0,constantes_new - constantes);	  
	  if(alpha .== .NaN){println(alpha);}

if(log(unif) .< alpha){		
	aceit[0][i] += 1;
	return dropr(X_new, (rows(X_new)-1));
}
else{							 		
	return dropr(X, (rows(X)-1)); 
}
	
	
}





cond_VR(const psi,const num_intervalos, const y, const W, X, const S,
		const B_matrix, const dist_inicial, aceit, total_testes,
		const intervalos_tempos, const obs_bloco,
		const k,const valores_r, const probs_r, const classe_priori_r, const hiperpar_r){
					
decl e;
e = columns(B_matrix);
							 				
X = X | (X[rows(X)-1][0] ~ S[1] ~ 0);		 	  

decl aux1, aux2, aux3, aux4, i;
		  
aux1 = round(range(0,(num_intervalos-1)).*((rows(W)-1) ./ num_intervalos)) ~ (rows(W)-1); 
aux3 = mincindex(cabs((W[aux1[0]][] - X[][1])') ');
				
for(i = 1; i < num_intervalos; ++i){
 
	aux3 |= mincindex(cabs((W[aux1[i]][] - X[][1])') ') .== aux3[i-1] .? mincindex(cabs((W[aux1[i]][] - X[][1])') ') + 1
																.: mincindex(cabs((W[aux1[i]][] - X[][1])') ');
}

aux3 |= rows(X)-1; 

aux4 = cumulate(countc(W, X[dropr(aux3,rows(aux3)-1)][1]));

decl X_interval, W_interval, V_b_cond, R_b_cond, delta, constantes, X_prop, constantes_prop, X_new, indices;

if(num_intervalos.== 1){	 

	delta = constant(S[0], 1, e);
	constantes = constantes_func(intervalos_tempos,obs_bloco, delta,
							psi,y,W,X[][0:1],X[][2],B_matrix,dist_inicial,
							e,k,valores_r,probs_r,classe_priori_r,hiperpar_r);	 																 				
	X_prop = gerando_vr(intervalos_tempos,obs_bloco,
								psi,y,W,B_matrix,dist_inicial,
				e,k,valores_r,probs_r,classe_priori_r,hiperpar_r);	
	constantes_prop = X_prop[1];
	X_prop = X_prop[0];
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ S[1] ~ 0);
	X_new = teste_mh(X, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], 0);	 
}
else{	

	i = 0; 
	
	X_interval = X[aux3[i]:aux3[i+1]][]; 	   
								  
	if(rows(X_interval) .== 1){
	X_new = X_interval;
	}
	else{
	W_interval = W[(aux4[i]-1):(aux4[i+1]-1)];   			
	
	V_b_cond = X_interval[rows(X_interval)-1][0]; 
	
	delta = constant(S[0], 1, e);

	decl aux5, intervalos_tempos_aux, obs_bloco_aux;
	aux5 = vecindex(intervalos_tempos[][0] .== W_interval[0]);
	intervalos_tempos_aux = intervalos_tempos[aux5:(aux5+rows(W_interval)-2)][];
	obs_bloco_aux = obs_bloco[aux5:(aux5+rows(W_interval)-1)][];
								
	constantes = constantes_func_cond(intervalos_tempos_aux, obs_bloco_aux, delta,
					psi,y, W_interval, X_interval[][0:1], X_interval[][2], B_matrix, dist_inicial,
					e, k, valores_r, probs_r, classe_priori_r,
					hiperpar_r);
	X_prop = gerando_vr_prim(intervalos_tempos_aux, obs_bloco_aux, V_b_cond,
						psi,y, W_interval, B_matrix, dist_inicial, e, k,
						valores_r, probs_r, classe_priori_r, hiperpar_r);		
	constantes_prop = X_prop[1];	  
	X_prop = X_prop[0];
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1:2]);
	X_new = teste_mh(X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);
	}
	
	for(i = 1; i < num_intervalos -1; ++i){   
	
		delta = constant(X[aux3[i]][1], 1, e);  
		delta[X_new[rows(X_new)-1][0]] = X_new[rows(X_new)-1][1]; 
		X_interval = X[aux3[i]:aux3[i+1]][];		
		if(rows(X_interval) .== 1){	 
			X_new = X_new | X_interval;
		}
		else{
		W_interval = W[(aux4[i]-1):(aux4[i+1]-1)];  
		V_b_cond = X_new[rows(X_new)-1][0] ~ X_interval[rows(X_interval)-1][0];

		decl aux5, intervalos_tempos_aux, obs_bloco_aux;
		aux5 = vecindex(intervalos_tempos[][0] .== W_interval[0]);
		intervalos_tempos_aux = intervalos_tempos[aux5:(aux5+rows(W_interval)-2)][];
		obs_bloco_aux = obs_bloco[aux5:(aux5+rows(W_interval)-1)][];
		obs_bloco_aux[0] = cumulate(obs_bloco)[aux5];
   																								
		constantes = constantes_func_cond(intervalos_tempos_aux, obs_bloco_aux, delta,
						psi,y, W_interval, X_interval[][0:1], X_interval[][2], B_matrix, B_matrix[V_b_cond[0]][],	 
						e, k, valores_r, probs_r, classe_priori_r,
						hiperpar_r);
		X_prop = gerando_vr_inter(intervalos_tempos_aux, obs_bloco_aux, V_b_cond, delta,
				psi, y, W_interval, B_matrix, dist_inicial, e, k,
				valores_r, probs_r, classe_priori_r, hiperpar_r);
														 
		constantes_prop = X_prop[1];
		X_prop = X_prop[0];		
		X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1:2]);   	
		X_new = X_new |teste_mh(X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);

		}

	}
	
	i = num_intervalos - 1;	
	
	delta = constant(X[aux3[i]][1], 1, e);
	delta[X_new[rows(X_new)-1][0]] = X_new[rows(X_new)-1][1];
	X_interval = X[aux3[i]:aux3[i+1]][];
	
	if(rows(X_interval) .== 1){
		X_new = X_new | X_interval;
	}
	else{				
	W_interval = W[(aux4[i]-1):(aux4[i+1]-1)]; 	 
	V_b_cond = X_new[rows(X_new)-1][0];
	R_b_cond = X_new[rows(X_new)-1][2];
	
	decl aux5, intervalos_tempos_aux, obs_bloco_aux;
	aux5 = vecindex(intervalos_tempos[][0] .== W_interval[0]);		
    intervalos_tempos_aux = intervalos_tempos[aux5:(rows(intervalos_tempos)-1)][];
	obs_bloco_aux = obs_bloco[aux5:(aux5+rows(W_interval)-1)][];	
	obs_bloco_aux[0] = cumulate(obs_bloco)[aux5];

	constantes = constantes_func(intervalos_tempos_aux, obs_bloco_aux, delta,
					psi,y, W_interval, X_interval[][0:1], X_interval[][2], B_matrix, B_matrix[V_b_cond[0]][],
					e, k, valores_r, probs_r, classe_priori_r,
					hiperpar_r);
	X_prop = gerando_vr_ultimo(intervalos_tempos_aux, obs_bloco_aux, V_b_cond, R_b_cond, delta,
			psi, y, W_interval, B_matrix, dist_inicial, e, k,
			valores_r, probs_r, classe_priori_r, hiperpar_r);
														 
	constantes_prop = X_prop[1];
	X_prop = X_prop[0];		
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1:2]);   	
	X_new = X_new |teste_mh(X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);

	}
	  
}		  	

if(X_new[rows(X_new)-1][1] .== 200){	
	X_new = X_new[0:(rows(X_new)-2)][];	
}

return X_new;

}




cond_theta_var(const k, const e, const X_aux, const gamma, const dir){

decl theta1, theta2_aux, theta2, j;
theta1 = constant(.NaN, 1, k);
theta2_aux = constant(.NaN, k, k);

decl aux3, mj, deltaj, Aj, Bj, aux4, aux1log, aux2log;

aux1log = constant(.NaN, 1, k); 
aux2log = constant(.NaN, 1, k);

for(j = 0; j < k; ++j){					 
				 
aux3 = X_aux[vecindex(X_aux[][0] .== j .|| X_aux[][0] .== j + k)][];  
mj = rows(aux3); 
deltaj = sumc(aux3)[1]; 

if(X_aux[rows(X_aux)-1][0] .== j .|| X_aux[rows(X_aux)-1][0] .== j + k){ 
Aj = gamma[j][0] + mj - 1; 
}
else{ 	  
Aj = gamma[j][0] + mj; 
}
Bj = gamma[j][1] + deltaj; 

theta1[j] = rangamma(1, 1, Aj, Bj);	

aux1log[j] = Aj - 1;  
aux2log[j] = Bj;  

}

decl logP2;
logP2 = sumr(log(theta1) .* aux1log) - sumr(theta1 .* aux2log);
theta1 = (theta1 ~ theta1)[0:(e-1)];

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
			
for(j = 0; j < k; ++j){

aux4 = randirichlet(1, (dir + mdots_aux)[j][]); 
theta2_aux[j][] = aux4 ~ (1-sumr(aux4));	
	
}	

theta2 = (setdiagonal(theta2_aux, 0) ~ diagonalize(theta2_aux)) |
		 (theta2_aux ~ constant(0, k, k));
theta2 = theta2[0:(e-1)][0:(e-1)];
theta2 = theta2 ./ sumr(theta2);   	 

decl logP3, aux3log;
aux3log = theta2;			 

if(rows(aux3log) != 2*k){ 	 
aux3log ~= constant(0,rows(aux3log),2*k-rows(aux3log));
}

aux3log = aux3log[][0:k-1] + aux3log[][k:2*k-1];
aux3log = aux3log[0:k-1][]; 		

aux3log = (mdots_aux + dir - 1) .* log(aux3log);	 

logP3 = 0;
for(j = 0; j < k; ++j){
logP3 += sumr(deletec(aux3log[j][]));
}

decl theta_Q;
theta_Q = setdiagonal(theta2 .* theta1', -theta1);	

return {theta_Q, logP2};

}


cond_theta_fixa(const k, const e, const X_aux, const Q){

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
								  
return {Q,logP2+logP3};

}



cond_theta(const k, const e, const X_aux, const Q, const gamma, const dir){

decl theta;
if(rows(gamma).==0) { 
	theta = cond_theta_fixa(k,e,X_aux,Q);
}
else{	
	theta = cond_theta_var(k,e,X_aux,gamma,dir);
}
	
return theta;

}



logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, j, R){	  

decl delta, integral, i, aux4, logvero;

delta = X[][1];
for(i = 1; i < e; ++i){
	delta ~= X[][1];
}

integral = int_g(psi,intervalos_tempos, delta, R); 
integral = selectrc(integral,range(0,rows(X)-1),X[][0])'; 	 	

decl soma_loglambda;
soma_loglambda = constant(.NaN, rows(X),1);

for(i = 0; i < rows(X); ++i){  

	if(obs_bloco[j[i]+1] .== 0){
		soma_loglambda[i][] = 0;
	}
	else{	
		aux4 = sumc(log(func_g(psi,y[cumulate(obs_bloco)[j[i]]:cumulate(obs_bloco)[j[i]+1]-1], delta[i][],R[i])));
		soma_loglambda[i][] = aux4[X[i][0]];		 
	}
}												  

logvero = sumc(-integral + soma_loglambda); 
			  
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




cond_Psi_MH(e, X,R, S, k, g_constante, psi, gamma_psi, matriz_cov, y, obs_bloco, intervalos_tempos, aceit, total_testes, prioris_psi, par_priori_psi){

decl i, j, ind, psi_new;	   

decl logP4 = 0;

// AMOSTRANDO DIRETO: CONJUGADA!!
for(j = 0; j < sizeof(g_constante); ++j){
	i = g_constante[j];  
	ind = vecindex(X[][0] .== i .|| X[][0] .== i + k); 
	if(sizeof(ind) .== 0){	 
		psi[i] = zeros(1,sizec(psi[i]));
	}
	else{				   
		decl deltai, ni;
		deltai = sumc(X[ind][3]);    
		ni = sumc(obs_bloco[ind+1][]);	

		psi[i] = rangamma(1, 1, gamma_psi[j][0] + ni, gamma_psi[j][1] + deltai);

		logP4 += -psi[i]*(gamma_psi[j][1]) + (gamma_psi[j][0] - 1).*log(psi[i]);
	}
}

// PASSO MH
decl dens_transicao, dens_transicao_new;		  

for(j = 0; j < sizeof(psi)-sizeof(g_constante); ++j){
			
	i = deletec(range(0,k-1), g_constante)[j]; 
	ind = vecindex(X[][0] .== i .|| X[][0] .== i + k);
	
	if(sizeof(ind) .== 0){	 
		psi[i] = zeros(1,sizec(psi[i]));
	}
	else{  			 
		psi_new = psi;								 
		psi_new[i] = rmvnorm(1, columns(psi[i]), psi[i], matriz_cov[j]);
									
		if(avaliar_lambda(X[][0:1] | (X[rows(X)-1][0] ~ S[1]), psi_new, e, R) > 0
					&& condicao_sinal_psi(psi_new) .== 1){																	
			dens_transicao = logverossimilhanca(e, psi,y, obs_bloco, X[ind][], intervalos_tempos[ind][], ind, R[ind])[0]  	  
									+ dens_priori_psi(psi, prioris_psi, par_priori_psi);
			dens_transicao_new = logverossimilhanca(e, psi_new,y, obs_bloco, X[ind][], intervalos_tempos[ind][], ind, R[ind])[0] 	
									+ dens_priori_psi(psi_new, prioris_psi, par_priori_psi);
			psi[i] = teste_mhPsi(psi[i], dens_transicao, psi_new[i], dens_transicao_new, &aceit[0], &total_testes[0], j);
			}
		else{	   			  
			total_testes[0][j] += 1;	
		}
	}
}

return {psi,logP4};

}





calib_MCMC(burnin, iterWVR, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, aceitCTMC,
				total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				aceitPsi, total_testesPsi, psi_aux2,valores_r, probs_r, classe_priori_r, hiperpar_r, prioris_psi, par_priori_psi){
					
decl psi_burnin;				
psi_burnin = constant(.NaN, burnin, sumr(psi_tamanho));

decl i, W, X_aux, Omega, intervalos_tempos, obs_bloco, psi_aux, j;

for(i = 0; i < burnin; ++i){

for(j = 0; j < iterWVR; ++j){
	W = cond_W(X[][0:1], S, Q_matrix);
	intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];					 
	obs_bloco = dropr(countc(y,W), rows(W));

	X = cond_VR(psi,num_intervalos,y,W, X, S,
		B_matrix, dist_inicial, &aceitCTMC[0], &total_testesCTMC[0],
		intervalos_tempos, obs_bloco,
		k, valores_r, probs_r, classe_priori_r, hiperpar_r);
}

	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]);
	Q_matrix = cond_theta(k, e, X_aux, Q_matrix, gamma, dir)[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];

	obs_bloco = countc(y,X[][1]);
	
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1],X[][2], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi[0], &total_testesPsi[0], prioris_psi, par_priori_psi)[0];

	psi_aux2 = <>;
	for(j = 0; j < sizeof(psi); ++j){
		psi[j] = psi_aux[j] == 0 .? psi[j] .: psi_aux[j];
		psi_aux2 ~= psi_aux[j];
	}
	
	psi_burnin[i][] = psi_aux2;	  
}	  
		  
return {X, Q_matrix, psi, psi_burnin};

}



calibracao(const_calib, const aceitPsi, const total_testesPsi){
						
decl k, prob_aceit;
					   
for(k = 0; k < columns(const_calib); ++k){

	prob_aceit = aceitPsi[k] ./ total_testesPsi[k];	 
	const_calib[0][k] = prob_aceit .< 0.2 .? const_calib[0][k] .* 0.9 .:
						prob_aceit .> 0.4 .? const_calib[0][k] .* 1.1 .: const_calib[0][k];
								
}

return const_calib;

}




gerar_processo(const S, const S_prev, const estados, const Q_matrix, const psi, const e,
					const valores_r, const probs_r, const grid, const final_X){

decl limite_trunc_exp;
limite_trunc_exp = {- final_X[][2] ./ psi[0][0], +.Inf};

decl Tj, Zj, Rj, j, X;

j = 0;

Tj = ranexpt(-Q_matrix[final_X[0]][final_X[0]],
					limite_trunc_exp[final_X[0]]-(S[1]-final_X[1]) );

if(columns(estados) .== 2){ 
   	Zj = dropc(estados, final_X[0]); 
}
else{
   	Zj = dropc(estados,final_X[0])[vecindex(ranmultinomial(1,dropc(Q_matrix[final_X[0]][],final_X[0]) ./ -Q_matrix[final_X[0]][final_X[0]])) .== 1];
}

Rj = valores_r[Zj][vecindex(ranmultinomial(1, probs_r[Zj]) .== 1)];   


while(Tj[j] .< S_prev){

	limite_trunc_exp = {- Rj[j] ./ psi[0][0], +.Inf};   
	Tj ~= Tj[j] + ranexpt(-Q_matrix[Zj[j]][Zj[j]], limite_trunc_exp[Zj[j]]);
	if(columns(estados) .== 2){ 
    	Zj ~= dropc(estados, Zj[j]); 
  	}
 	else{
   		Zj ~= dropc(estados,Zj[j])[vecindex(ranmultinomial(1,dropc(Q_matrix[Zj[j]][],Zj[j]) ./ -Q_matrix[Zj[j]][Zj[j]])) .== 1];
  	}
	Rj ~= valores_r[Zj[j+1]][vecindex(ranmultinomial(1, probs_r[Zj[j+1]]) .== 1)];
	j += 1;

}

X = Zj' ~ Tj' ~ Rj';	

X[][1] += S[1];

X = final_X | X;
X = X[vecindex((X[][1] .<= S[1]+S_prev) .== 1)][];   	

decl lambda_max;
lambda_max = max( avaliar_lambda(X[][0:1] | (X[rows(X)-1][0] ~ (S[1]+S_prev)), psi, e, X[][2]));

decl N, PP, prob, u, NHPP;

N = ranpoisson(1,1,(S_prev+S[1]-X[0][1])*lambda_max); 	  
PP = sortr(X[0][1]+ranu(1,N)*(S_prev+S[1]-X[0][1]));  

decl aux1, aux3, aux4, aux, lambda;
aux1 = sumc(X[][1] .<= PP)-1;										  
aux3 = <>;
aux4 = <>;
for(j = 0;j<e;++j){
	aux3 ~= X[aux1][1];
	aux4 ~= X[aux1][2];
}					  
// gerando eventos NHPP			
aux = func_g(psi, PP',aux3,aux4); 
prob = selectrc(aux,range(0,columns(PP)-1, 1)', X[aux1][0]') ./ lambda_max;
u = ranu(1,N);
NHPP = PP[vecindex(u .<= prob)];
NHPP = dropc(NHPP,vecindex(NHPP .<= S[1]));    


// guardando grid da FI
aux1 = sumc(X[][1] .<= grid)-1;									  
aux3 = <>;
aux4 = <>;
for(j = 0;j<e;++j){
	aux3 ~= X[aux1][1];
	aux4 ~= X[aux1][2];
}
lambda = func_g(psi, grid',aux3,aux4);
lambda = selectrc(lambda,range(0,columns(grid)-1, 1)', X[aux1][0]');	


decl delta, integral, integral2, i;
delta = X[][1];
for(i = 1; i < e; ++i){
	delta ~= X[][1];
}

integral = int_g(psi,(X[][1]|(S[1]+S_prev))[:(rows(X)-1)] ~ (X[][1]|(S[1]+S_prev))[1:], delta,X[][2]); 
integral = sumc(selectrc(integral,range(0,rows(X)-1),X[][0])');

integral2 = int_g(psi, X[0][1] ~ S[1], delta[0][],X[0][2]);	  
integral2 = integral2[X[0][0]];	
integral = integral - integral2;


return{lambda,columns(NHPP), integral};	  

}
