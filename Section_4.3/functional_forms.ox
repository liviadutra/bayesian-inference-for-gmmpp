const decl
c1 = 1,
c2 = 1;


func_g(psi,s,delta,R){
	
	decl g1, g2, g3;
	g1 = psi[0][0] + psi[0][1] .* (s - delta[][0]) .^ c1; 	
	g2 = R[][0] + psi[1][0] .* (s - delta[][1]) .^ c2;
	return g1 ~ g2;
}

int_g(psi,intervalo, delta,R){

	decl s1, s2, int1, int2, int3, int4;
	s1 = intervalo[][0];
	s2 = intervalo[][1];

	int1 = psi[0][0] .* (s2 - s1) + (psi[0][1]/(c1+1)) .*
				((s2 - delta[][0]) .^ (c1+1) - (s1 - delta[][0]) .^ (c1+1));
	int2 = R[][0] .* (s2 - s1) + (psi[1][0]/(c2+1)) .*
				((s2 - delta[][1]) .^ (c2+1) - (s1 - delta[][1]) .^ (c2+1));
	return int1 ~ int2;
}

condicao_sinal_psi(psi){
decl teste;
teste = psi[0][1] .< 0 && psi[1][0] .> 0;

return teste;

}