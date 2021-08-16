#include <oxprob.h>
#include <oxfloat.oxh>

func_g(psi,s,T,b2){
	decl g1, gamma_aux, g2;
	g1 = psi[0][1] + psi[0][0].*
				probn( psi[0][3]+sqrt(2 .* M_PI)./psi[0][0].*psi[0][2].*s ); 
	gamma_aux = (psi[0][1] + psi[0][0]*probn(psi[0][3] + sqrt(2*M_PI)/psi[0][0] * psi[0][2] * T) - b2)/
						(probn(psi[1][1]));
	g2 = b2 + gamma_aux*(probn(psi[1][1] - sqrt(2*M_PI)/gamma_aux *psi[1][0] * (s-T)) );
	return g1 ~ g2;
}

int_g(psi,T,S,b2){
	decl int1, int2, gamma_aux, grid_int, tam_grid_int;
	tam_grid_int = 100;
	grid_int = range(0,T,T/tam_grid_int);
	int1 = func_g(psi,meanr(grid_int[0:(tam_grid_int - 1)]' ~ grid_int[1:tam_grid_int]'),T,b2)[][0];
	int1 = sumc(int1 .* (T/tam_grid_int));

	gamma_aux = (psi[0][1] + psi[0][0]*probn(psi[0][3] + sqrt(2*M_PI)/psi[0][0] * psi[0][2] * T) - b2)./(probn(psi[1][1]));
	grid_int = range(T,S[1],(S[1]-T)/tam_grid_int);
	int2 = 	func_g(psi,meanr(grid_int[0:(tam_grid_int - 1)]' ~ grid_int[1:tam_grid_int]'),T,b2)[][1];
	int2 = sumc(int2 .* ((S[1]-T)/tam_grid_int));
			 
	return int1 ~ int2;
}


condicao_sinal_psi(psi){
decl teste;
teste = psi[0][0] .> 0 .&& psi[0][1] .>= 0 .&& psi[0][2] .> 0 .&&
	 	 psi[1][0] .> 0	
		 .&& psi[1][1] .> 0 .&& psi[1][1] .< 3;

return teste;

}