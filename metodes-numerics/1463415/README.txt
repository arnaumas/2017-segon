Per a compilar cada programa cal executar 'make <Nom del programa>. 'make all' compila tots els programes alhora.

El programa 'integralExp' calcula la integral de la funció exp(-x^2)/sqrt(1 - x^2) entre -1 i 1 mitjançant Gauss-Legendre i Gauss-Chebyshev. Cal donar el nombre de nodes que es faran servir. Si per exemple volem calcular les integrals amb 8 nodes:
	> make integralExp
	> ./integralExp
	>> Introduïu el nombre de nodes: 8
	>> Quadratura per Gauss-Chebyshev = 2.026438
	>> Quadratura per Gauss-Legendre = 1.951603

El programa 'integralAbs' calcula la integral de la funció |x| entre -1 i 1 mitjançant Gauss-Legendre, Gauss-Chebyshev i per la regla dels trapezis. Cal donar el nombre de nodes que es faran servir. Si per exemple volem calcular les integrals amb 8 nodes:
	> make integralAbs
	> ./integralAbs
	>> Introduïu el nombre de nodes: 8
	>> Quadratura per Gauss-Chebyshev = 1.026172
	>> Quadratura per Gauss-Legendre = 1.011528
	>> Quadratura per Trapezis = 1.020408

El programa 'integralElipse' calcula la longitud d'arc de l'el·lipse d'equació (x/2)^2 + (2y)^2 = 1 entre -1 i 1.mitjançant Gauss-Legendre, Gauss-Chebyshev i per la regla dels trapezis. Cal donar el nombre de nodes que es faran servir. Si per exemple volem calcular les integrals amb 8 nodes:
	> make integralElipse
	> ./integralElipse
	>> Introduïu el nombre de nodes: 8
	>> Quadratura per Gauss-Chebyshev = 2.019193
	>> Quadratura per Gauss-Legendre = 2.006145
	>> Quadratura per Trapezis = 2.006516

