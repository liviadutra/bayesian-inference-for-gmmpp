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



ranexpt(const lambda, const b){

	decl x;
	x = -(1 ./ lambda) .* log(1 - ranu(1,1).*(1 - exp(-lambda .* b)) );
  
	return x;			 
}



avaliar_lambda(X, psi, e){	  

	decl aux, epsilon, delta, i, inicial, final;
	aux = rows(X)-2;
	epsilon = 10^(-10); 

	delta = X[][1];
	for(i = 0; i < e; ++i){
		delta ~= X[][1];
	}
	delta = delta[:aux][]; 
	
	inicial = func_g(psi, X[:aux][1], delta);	 
	inicial = selectrc(inicial,range(0,aux),X[:aux][0])'; 
	
	final = func_g(psi, X[1:][1]-epsilon, delta);	
	final = selectrc(final,range(0,aux),X[:aux][0])'; 	 

	return inicial ~ final;		 
}

  

sdtmc(psi,const y, const W, const B_matrix, const dist_inicial, const e){	 

decl  intervalos_tempos, obs_bloco, V, integral, avaliar_func, soma_loglambda, i, delta, j, constantes_new;

if(rows(W) .== 2){
	intervalos_tempos =	W[0] ~ W[1];
}
else{
	intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}
obs_bloco = dropr(countc(y,W), rows(W));
V = constant(.NaN, rows(W)-1, 1);

// para amostrar o primeiro estado
i = 0;
delta = constant(W[0], 1, e);

avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);	 
integral = int_g(psi,intervalos_tempos[i][], delta);

if(obs_bloco[i+1] .== 0){
	soma_loglambda = constant(0, 1, e);	 
}
else{	
	soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));	 
}	  

decl aux3, aux4;
aux3 = -integral + soma_loglambda + log(dist_inicial);	
aux3[vecindex(avaliar_func)] = -.Inf;											
aux4 = maxr(aux3);		
constantes_new = 0;
V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);

decl delta_aux;
delta_aux = delta[V[i]];
delta = constant(W[i+1], 1, e);
delta[V[i]] = delta_aux;

// para amostrar os outros estados									
for(i = 1; i < rows(W)-1; ++i){

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0); 
	integral = int_g(psi,intervalos_tempos[i][], delta);

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]);
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4;			  		
	V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);
											
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;

}	

V = (V | V[rows(V)-1]) ~ W;	  

decl indices, X;
indices = 1 | (V[1:(rows(V)-1)][0] .== V[0:(rows(V)-2)][0] .? 0 .: 1);	
X = V[vecindex(indices .== 1)][];	
	 
return{X, constantes_new};

}


sdtmc_prim(psi,const y, const W, const V_b_cond, const B_matrix, const dist_inicial, const e){
				  
decl  intervalos_tempos, obs_bloco, V, integral, avaliar_func, soma_loglambda, i, delta, constantes_new;

if(rows(W) .== 2){
intervalos_tempos =	W[0] ~ W[1];
}
else{
intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}
obs_bloco = dropr(countc(y,W), rows(W));
V = constant(.NaN, rows(W), 1);
V[rows(W)-1][0] = V_b_cond;

i = 0;
delta = constant(W[0], 1, e); 
avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
integral = int_g(psi,intervalos_tempos[i][], delta);

if(obs_bloco[i+1] .== 0){
	soma_loglambda = constant(0, 1, e);		
}
else{	
	soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));		 
}	  

decl aux3, aux4;
aux3 = -integral + soma_loglambda + log(dist_inicial);		  
aux3[vecindex(avaliar_func)] = -.Inf;
aux4 = maxr(aux3);
constantes_new = 0;	 
V[i] = vecindex(ranmultinomial(1,  exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);  
					 			   		
decl delta_aux;
delta_aux = delta[V[i]];
delta = constant(W[i+1], 1, e);	
delta[V[i]] = delta_aux;	  
								
for(i = 1; i < rows(W)-2; ++i){	

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);	   

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]);
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4;
	V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);
										
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;					

}  

i = rows(W)-2;

if(i .== 0){}
else{	 
	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1][0]][]) + log(B_matrix[][V[i+1][0]]');
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);			   
		
	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4; 
}

decl prob;				
prob = exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4);	 
prob[V[rows(W)-1]] = 0;	
V[i] = vecindex(ranmultinomial(1, prob ./ sumr(prob)) .== 1); 
V = V ~ W;
						 				
decl indices, X;
indices = 1 | (V[1:(rows(V)-1)][0] .== V[0:(rows(V)-2)][0] .? 0 .: 1);	
X = V[vecindex(indices .== 1)][];
								
return {dropr(X, rows(X)-1), constantes_new};	 

}


sdtmc_inter(psi,const y, const W, const V_b_cond, const B_matrix, const e, delta){	  
					 							  
decl  intervalos_tempos, obs_bloco, V, i, integral, avaliar_func, soma_loglambda, constantes_new;

if(rows(W) .== 2){
intervalos_tempos =	W[0] ~ W[1];
}
else{
intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}

obs_bloco = dropr(countc(y,W), rows(W));
											  
V = constant(.NaN, rows(W), 1); 
V[rows(W)-1][0] = V_b_cond[1];   

i = 0;

avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
integral = int_g(psi,intervalos_tempos[i][], delta);

if(obs_bloco[i+1] .== 0){
	soma_loglambda = constant(0, 1, e);
}
else{	
	soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta))); 
}	  

decl aux3, aux4;
aux3 = -integral + soma_loglambda + log(B_matrix[V_b_cond[0]][]);		  
aux3[vecindex(avaliar_func)] = -.Inf;
aux4 = maxr(aux3);

constantes_new = 0;	
V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);
								
decl delta_aux;			
delta_aux = delta[V[i]];  
delta = constant(W[i+1], 1, e);	 
delta[V[i]] = delta_aux;	  	

for(i = 1; i < rows(W)-2; ++i){	

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);	   

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]);
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);

	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4;
	V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);
										
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;					
						
}  

i = rows(W)-2;	
if(i .== 0){}
else{	 

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0); 	
	integral = int_g(psi,intervalos_tempos[i][], delta);

	if(obs_bloco[i+1] .== 0){			  
		soma_loglambda = constant(0, 1, e);
	}
	else{		  
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1][0]][]) + log(B_matrix[][V[i+1][0]]');
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);

	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4;
  					
}				  

decl prob;
prob = exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4);		 
prob[V[rows(W)-1]] = 0; 	 					
		  
V[i] = vecindex(ranmultinomial(1, prob ./ sumr(prob)) .== 1);	
V = V ~ W;	 	 	  				
						
decl indices, X;
indices = 1 | (V[1:(rows(V)-1)][0] .== V[0:(rows(V)-2)][0] .? 0 .: 1);	
X = V[vecindex(indices .== 1)][]; 

return {dropr(X, rows(X)-1), constantes_new};

}



sdtmc_ultimo(psi,const y, const W, const V_b_cond, const B_matrix, const e, delta){		 
									
decl  intervalos_tempos, obs_bloco, V, i, integral, avaliar_func, soma_loglambda, constantes_new;

if(rows(W) .== 2){
intervalos_tempos =	W[0] ~ W[1];
}
else{
intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}
obs_bloco = dropr(countc(y,W), rows(W));	
V = constant(.NaN, rows(W)-1, 1);

i = 0;
avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
integral = int_g(psi,intervalos_tempos[i][], delta);	 
if(obs_bloco[i+1] .== 0){
	soma_loglambda = constant(0, 1, e);	 	
}
else{	
	soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));		 
}	  

decl aux3, aux4;
aux3 = -integral + soma_loglambda + log(B_matrix[V_b_cond][]);		  
aux3[vecindex(avaliar_func)] = -.Inf;
aux4 = maxr(aux3);

constantes_new = 0;	 
V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);  
				 
decl delta_aux;
delta_aux = delta[V[i]];
delta = constant(W[i+1], 1, e);	
delta[V[i]] = delta_aux;	  

for(i = 1; i < rows(W)-1; ++i){	 				

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);	   
  
	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]);
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);		
	
	constantes_new += log(sumr(exp(aux3 - aux4))) + aux4;					   
	V[i] = vecindex(ranmultinomial(1, exp(aux3 - log(sumr(exp(aux3 - aux4))) - aux4)) .== 1);
							
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;					

}				
  														
return{V ~ dropr(W, (rows(W)-1)), constantes_new};	// eliminando o tempo final.

}



constantes_func_cond(psi,const y, const W, const X, const B_matrix, delta){	 

decl  e, intervalos_tempos, obs_bloco, V, integral, avaliar_func, soma_loglambda,
					i, constantes;

e = columns(B_matrix);
 
if(rows(W) .== 2){	
intervalos_tempos =	W[0] ~ W[1];
}
else{			if(rows(W).==1){println(W);}
intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}
obs_bloco = dropr(countc(y,W), rows(W));		
V = constant(.NaN, rows(W)-1, 1);

i = 0;
constantes = 0;	 
V[i] = X[i][0];	 

decl delta_aux;
delta_aux = delta[V[i]];		 		  			
delta = constant(W[i+1], 1, e);					
delta[V[i]] = delta_aux;				  				

decl aux3, aux4;									
for(i = 1; i < rows(W)-2; ++i){	

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);	

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]); 	 
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
		
	constantes += log(sumr(exp(aux3 - aux4))) + aux4;							
	V[i] = X[maxc(vecindex(W[i] .>= X[][1]))][0];  	 
											
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;			  

}	  					

i = rows(W)-2;	
if(i .== 0){ }		
else{

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);

	if(obs_bloco[i+1] .== 0){			  
		soma_loglambda = constant(0, 1, e);
	}
	else{		  
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}
								 
	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1][0]][]) + log(B_matrix[][X[rows(X)-1][0]]');
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
	
	constantes += log(sumr(exp(aux3 - aux4))) + aux4;
}

return constantes;

}



constantes_func(psi,const y, const W, const X, const B_matrix, delta){	   

decl  e, intervalos_tempos, obs_bloco, V, integral, avaliar_func, soma_loglambda,
					i, constantes;

e = columns(B_matrix);
 
if(rows(W) .== 2){
intervalos_tempos =	W[0] ~ W[1];
}
else{
intervalos_tempos = W[:(rows(W)-2)] ~ W[1:];
}
obs_bloco = dropr(countc(y,W), rows(W));
V = constant(.NaN, rows(W)-1, 1);   			 				 

i = 0;
constantes = 0;	
V[i] = X[i][0];

decl delta_aux;
delta_aux = delta[V[i]];		 		  			
delta = constant(W[i+1], 1, e);					
delta[V[i]] = delta_aux;			  				

decl aux3, aux4;					
for(i = 1; i < rows(W)-2; ++i){		

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);

	if(obs_bloco[i+1] .== 0){
		soma_loglambda = constant(0, 1, e);
	}
	else{
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}						   
						 
	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1]][]); 
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
	
	constantes += log(sumr(exp(aux3 - aux4))) + aux4;
	V[i] = X[maxc(vecindex(W[i] .>= X[][1]))][0];  	 
											
	delta_aux = delta[V[i]];
	delta = constant(W[i+1], 1, e);
	delta[V[i]] = delta_aux;	 

}	  					

i = rows(W)-2;		
if(i .== 0){}	
else{

	avaliar_func = (func_g(psi, W[i], delta) .< 0) .|| (func_g(psi, W[i+1], delta) .< 0);
	integral = int_g(psi,intervalos_tempos[i][], delta);  

	if(obs_bloco[i+1] .== 0){			  
		soma_loglambda = constant(0, 1, e);
	}
	else{		  
		soma_loglambda = sumc(log(func_g(psi,y[cumulate(obs_bloco)[i]:cumulate(obs_bloco)[i+1]-1], delta)));
	}

	aux3 = -integral + soma_loglambda + log(B_matrix[V[i-1][0]][]);	 
	aux3[vecindex(avaliar_func)] = -.Inf;
	aux4 = maxr(aux3);
		
	constantes += log(sumr(exp(aux3 - aux4))) + aux4; // a constante é calculada antes de obter a condição da indicadora!!
				  	 		
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



teste_mh(const psi, const X, const constantes, const X_new, const constantes_new, const aceit, const total_testes, const i){
				   			 
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




cond_V(const psi,const num_intervalos, const y, const W, X, const S, const B_matrix, const dist_inicial, const aceit, const total_testes){
				
decl e;
e = columns(B_matrix);

X = X | (X[rows(X)-1][0] ~ S[1]);	 	  

decl aux1, aux3, aux4, i;

aux1 = round(range(0,(num_intervalos-1)).*((rows(W)-1) ./ num_intervalos)) ~ (rows(W)-1); 
aux3 = mincindex(cabs((W[aux1[0]][] - X[][1])') ');
				
for(i = 1; i < num_intervalos; ++i){
 
	aux3 |= mincindex(cabs((W[aux1[i]][] - X[][1])') ') .== aux3[i-1] .? mincindex(cabs((W[aux1[i]][] - X[][1])') ') + 1
																.: mincindex(cabs((W[aux1[i]][] - X[][1])') ');
}

aux3 |= rows(X)-1; 

aux4 = cumulate(countc(W, X[dropr(aux3,rows(aux3)-1)][1]));

decl X_interval, W_interval, V_b_cond, delta, constantes, X_prop, constantes_prop, X_new, indices;

if(num_intervalos.== 1){	 

	delta = constant(S[0], 1, e);
	constantes = constantes_func(psi,y, W, X, B_matrix, delta);	
																 			
	X_prop = sdtmc(psi,y, W, B_matrix, dist_inicial, e);	
	constantes_prop = X_prop[1];
	X_prop = X_prop[0];
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ S[1]); 
															
	X_new = teste_mh(psi, X, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], 0);	 
		
}				 
else{		 

	i = 0; 
	
	X_interval = X[aux3[i]:aux3[i+1]][];  	   
	W_interval = W[(aux4[i]-1):(aux4[i+1]-1)];  
	V_b_cond = X_interval[rows(X_interval)-1][0]; 
	delta = constant(S[0], 1, e);
	
	constantes = constantes_func_cond(psi,y, W_interval, X_interval, B_matrix, delta);
	X_prop = sdtmc_prim(psi,y, W_interval, V_b_cond, B_matrix, dist_inicial, e);		
	constantes_prop = X_prop[1];	  
	X_prop = X_prop[0];
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1]);

	X_new = teste_mh(psi, X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);
	
	for(i = 1; i < num_intervalos -1; ++i){ 
	
		delta = constant(X[aux3[i]][1], 1, e); 
		delta[X_new[rows(X_new)-1][0]] = X_new[rows(X_new)-1][1]; 
				
		X_interval = X[aux3[i]:aux3[i+1]][];
		
		W_interval = W[(aux4[i]-1):(aux4[i+1]-1)];  			
											
		V_b_cond = X_new[rows(X_new)-1][0] ~ X_interval[rows(X_interval)-1][0];	
		constantes = constantes_func_cond(psi,y, W_interval, X_interval, B_matrix, delta);
		
		X_prop = sdtmc_inter(psi,y, W_interval, V_b_cond, B_matrix, e, delta); 
										
		constantes_prop = X_prop[1];	  
		X_prop = X_prop[0];			  
		X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1]);

		X_new = X_new | teste_mh(psi, X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);	  												 

		indices = 1 | (X_new[1:(rows(X_new)-1)][0] .== X_new[0:(rows(X_new)-2)][0] .? 0 .: 1);	
		X_new = X_new[vecindex(indices .== 1)][];
	
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
	
	constantes = constantes_func(psi,y, W_interval, X_interval, B_matrix, delta);
																 			
	X_prop = sdtmc_ultimo(psi,y, W_interval, V_b_cond, B_matrix, e, delta);
	
	constantes_prop = X_prop[1];
	X_prop = X_prop[0];		
	X_prop = X_prop | (X_prop[rows(X_prop)-1][0] ~ X_interval[rows(X_interval)-1][1]);
	
	X_new = X_new | teste_mh(psi, X_interval, constantes, X_prop, constantes_prop, &aceit[0], &total_testes[0], i);
	}
	
	indices = 1 | (X_new[1:(rows(X_new)-1)][0] .== X_new[0:(rows(X_new)-2)][0] .? 0 .: 1);
	X_new = X_new[vecindex(indices .== 1)][];
		  
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
								  
return {theta_Q, logP2+logP3};

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



logverossimilhanca(e, psi,y, obs_bloco, X, intervalos_tempos, j){	  

decl delta, integral, i, aux4, logvero;

delta = X[][1];
for(i = 0; i < e; ++i){
	delta ~= X[][1];
}

integral = int_g(psi,intervalos_tempos, delta); 
integral = selectrc(integral,range(0,rows(X)-1),X[][0])'; 	 	

decl soma_loglambda;
soma_loglambda = constant(.NaN, rows(X),1);
												
for(i = 0; i < rows(X); ++i){ 

	if(obs_bloco[j[i]+1] .== 0){ 
		soma_loglambda[i][] = 0;
	}
	else{	
		aux4 = sumc(log(func_g(psi,y[cumulate(obs_bloco)[j[i]]:cumulate(obs_bloco)[j[i]+1]-1], delta[i][])));
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



cond_Psi_MH(e, X, S, k, g_constante, psi, gamma_psi, matriz_cov, y, obs_bloco, intervalos_tempos, aceit, total_testes, prioris_psi, par_priori_psi){

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
		deltai = sumc(X[ind][2]);    
		ni = sumc(obs_bloco[ind+1][]);	

		psi[i] = rangamma(1, 1, gamma_psi[j][0] + ni, gamma_psi[j][1] + deltai);

		logP4 += -psi[i]*(gamma_psi[j][1]) + (gamma_psi[j][0] - 1).*log(psi[i]);
	}
}

// PASSO MH
decl dens_transicao, dens_transicao_new;		  

for(j = 0; j < k-sizeof(g_constante); ++j){

	i = deletec(range(0,k-1), g_constante)[j]; 
	ind = vecindex(X[][0] .== i .|| X[][0] .== i + k); 	

	if(sizeof(ind) .== 0){
		psi[i] = zeros(1,sizec(psi[i]));
	}
	else{  			 
		psi_new = psi;									 
		psi_new[i] = rmvnorm(1, columns(psi[i]), psi[i], matriz_cov[j]);
									
		if(avaliar_lambda(X[][0:1] | (X[rows(X)-1][0] ~ S[1]), psi_new, e) > 0 &&
					condicao_sinal_psi(psi_new) .== 1){																	
			dens_transicao = logverossimilhanca(e, psi,y, obs_bloco, X[ind][], intervalos_tempos[ind][], ind)[0] 	
							 + dens_priori_psi(psi, prioris_psi, par_priori_psi);
			dens_transicao_new = logverossimilhanca(e, psi_new,y, obs_bloco, X[ind][], intervalos_tempos[ind][], ind)[0]	
					   		 + dens_priori_psi(psi_new, prioris_psi, par_priori_psi);
			psi[i] = teste_mhPsi(psi[i], dens_transicao, psi_new[i], dens_transicao_new, &aceit[0], &total_testes[0], j);
		}
		else{	   			  
			total_testes[0][j] += 1;	
		}
	}
}

logP4 += dens_priori_psi(psi, prioris_psi, par_priori_psi);

return {psi,logP4};

}


calib_MCMC(calib, iterWV, X, S, Q_matrix, psi, num_intervalos, y, B_matrix, dist_inicial, aceitCTMC,
				total_testesCTMC, k, g_constante, e, gamma, dir, gamma_psi, psi_tamanho, matriz_cov,
				aceitPsi, total_testesPsi, psi_aux2, prioris_psi, par_priori_psi){

decl psi_burnin;				
psi_burnin = constant(.NaN, calib, sumr(psi_tamanho));

decl i, W, X_aux, Omega, intervalos_tempos, obs_bloco, psi_aux, j;

for(i = 0; i < calib; ++i){	 

	for(j = 0; j < iterWV; ++j){
		W = cond_W(X, S, Q_matrix);   
		X = cond_V(psi, num_intervalos, y, W, X, S, B_matrix, dist_inicial, &aceitCTMC[0], &total_testesCTMC[0]);
	}
	X_aux = X[][0] ~ ((X[][1] | S[1])[1:] - X[:(rows(X)-1)][1]); 

	Q_matrix = cond_theta(k, e, X_aux, Q_matrix, gamma, dir)[0];
	Omega = 2*diagonal(-Q_matrix)';			  
	B_matrix = unit(e,e) + (1 ./ Omega) .* Q_matrix;

	intervalos_tempos = (X[][1]|S[1])[:(rows(X)-1)] ~ (X[][1]|S[1])[1:];
	obs_bloco = countc(y,X[][1]);
	
	psi_aux = cond_Psi_MH(e, X ~ X_aux[][1], S, k, g_constante, psi, gamma_psi, matriz_cov,
							y, obs_bloco, intervalos_tempos, &aceitPsi[0], &total_testesPsi[0],
							prioris_psi, par_priori_psi)[0];
		
	psi_aux2 = <>;
	for(j = 0; j < k; ++j){
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
					const grid, const final_X, const limite_trunc_exp){
				
// gerando uma CTMC X
decl Tj, Zj, j, X;

j = 0;
Tj = ranexpt(-Q_matrix[final_X[0]][final_X[0]],
					limite_trunc_exp[final_X[0]]-(S[1]-final_X[1]) );

if(columns(estados) .== 2){ 
   	Zj = dropc(estados, final_X[0]); 
}
else{
   	Zj = dropc(estados,final_X[0])[vecindex(ranmultinomial(1,dropc(Q_matrix[final_X[0]][],final_X[0]) ./ -Q_matrix[final_X[0]][final_X[0]])) .== 1];
}
while(Tj[j] .< S_prev){
	Tj ~= Tj[j] + ranexpt(-Q_matrix[Zj[j]][Zj[j]], limite_trunc_exp[Zj[j]]);
	if(columns(estados) .== 2){ 
    	Zj ~= dropc(estados, Zj[j]); 
  	}
 	else{
   		Zj ~= dropc(estados,Zj[j])[vecindex(ranmultinomial(1,dropc(Q_matrix[Zj[j]][],Zj[j]) ./ -Q_matrix[Zj[j]][Zj[j]])) .== 1];
  	}
	j += 1;
}

X = Zj' ~ Tj';	 
X[][1] += S[1];		
X = final_X | X;	
X = X[vecindex((X[][1] .<= S[1]+S_prev) .== 1)][];   

decl lambda_max;
lambda_max = max(avaliar_lambda(X[][0:1] | (X[rows(X)-1][0] ~ (S[1]+S_prev)), psi, e));
										
decl N, PP, prob, u, NHPP;

N = ranpoisson(1,1,(S_prev+S[1]-X[0][1])*lambda_max); 	  
PP = sortr(X[0][1]+ranu(1,N)*(S_prev+S[1]-X[0][1])); 

decl aux1, aux3, aux4, aux, lambda;
aux1 = sumc(X[][1] .<= PP)-1;										  
aux3 = <>;
aux4 = <>;
for(j = 0;j<e;++j){
	aux3 ~= X[aux1][1];
}					  
// gerando eventos NHPP			
aux = func_g(psi, PP',aux3); 
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
}
lambda = func_g(psi, grid',aux3);
lambda = selectrc(lambda,range(0,columns(grid)-1, 1)', X[aux1][0]');

decl delta, integral, integral2, i;
delta = X[][1];
for(i = 1; i < e; ++i){
	delta ~= X[][1];
}					

integral = int_g(psi,(X[][1]|(S[1]+S_prev))[:(rows(X)-1)] ~ (X[][1]|(S[1]+S_prev))[1:], delta); 
integral = sumc(selectrc(integral,range(0,rows(X)-1),X[][0])');

integral2 = int_g(psi, X[0][1] ~ S[1], delta[0][]);	  
integral2 = integral2[X[0][0]];	
integral = integral - integral2;

return{lambda,columns(NHPP), integral};

} 
