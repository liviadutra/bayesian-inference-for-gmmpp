const decl
c1 = 1,	
c2 = 1,
b3 = 0,
c3 = 1;


func_g(psi,s,delta){
	decl g1, g2, g3, g4, g5;
	g1 = psi[0][0] + psi[0][1] .* (s - delta[][0]) .^ c1; 
	g2 = psi[1][0] + psi[1][1] .* (s - delta[][1]) .^ c2;
	g3 = psi[2][0] + b3 .* (s - delta[][2]) .^ c3; 
	return g1 ~ g2 ~ g3;
}

int_g(psi,intervalo, delta){
	decl s1, s2, int1, int2, int3, int4, int5;
	s1 = intervalo[][0];
	s2 = intervalo[][1];
	int1 = psi[0][0] .* (s2 - s1) + (psi[0][1]/(c1+1)) .*
				((s2 - delta[][0]) .^ (c1+1) - (s1 - delta[][0]) .^ (c1+1));
	int2 = psi[1][0] .* (s2 - s1) + (psi[1][1]/(c2+1)) .*
				((s2 - delta[][1]) .^ (c2+1) - (s1 - delta[][1]) .^ (c2+1));
	int3 = psi[2][0] .* (s2 - s1) + (b3/(c3+1)) .*
				((s2 - delta[][2]) .^ (c3+1) - (s1 - delta[][2]) .^ (c3+1));
	return int1 ~ int2 ~ int3;
}


condicao_sinal_psi(psi){
decl teste;
teste = psi[0][0] .> 0 && psi[0][1] .> 0 &&
		psi[1][0] .> 0 && psi[1][1] .< 0;

return teste;

}
