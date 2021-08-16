const decl
//a1 = 1,
//b1 = 0.2,
c1 = 1,
c2 = 1,
b3 = 0,
c3 = 1;


func_g(psi,s,delta,R){
	
	decl g1, g2, g3;
	g1 = R[][0] + psi[0][0] .* (s - delta[][0]) .^ c1; 	
	g2 = R[][0] + psi[1][0] .* (s - delta[][1]) .^ c2;		
	g3 = R[][0] + b3 .* (s - delta[][2]) .^ c3;
	return g1 ~ g2 ~ g3;
}

int_g(psi,intervalo, delta,R){

	decl s1, s2, int1, int2, int3, int4;
	s1 = intervalo[][0];
	s2 = intervalo[][1];

	int1 = R[][0] .* (s2 - s1) + (psi[0][0]/(c1+1)) .*
				((s2 - delta[][0]) .^ (c1+1) - (s1 - delta[][0]) .^ (c1+1));
	int2 = R[][0] .* (s2 - s1) + (psi[1][0]/(c2+1)) .*
				((s2 - delta[][1]) .^ (c2+1) - (s1 - delta[][1]) .^ (c2+1));
	int3 = R[][0] .* (s2 - s1) + (b3/(c3+1)) .*
				((s2 - delta[][2]) .^ (c3+1) - (s1 - delta[][2]) .^ (c3+1));
	return int1 ~ int2 ~ int3;
}

condicao_sinal_psi(psi){
decl teste;
teste = psi[0][0] .< 0 && psi[1][0] .> 0;

return teste;

}




/*

função da reta:
f(s, delta) = a + b*(s - delta)^c

função exponencial:
f(s, delta) = a + b*exp{c*(s - delta)}



integral entre (s1, s2) de uma função do tipo: a + b*(s - delta)^c:
a*(s2 - s1) + (b/(c+1))*((s2 - delta)^(c+1) - (s1 - delta)^(c+1));

integral entre (s1, s2) de uma função do tipo: a + b*exp{c*(s - delta)}:
a*(s2 - s1) + (b/c)*exp{-delta*c}*(exp{c*s2} - exp{c*s1})

*/
