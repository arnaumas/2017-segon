#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"polinomis.h"

int main() {
	// Fem el càlcul de la integral de abs(x) entre -1 i 1 amb Gauss-Chebyshev
	int n;
	double* cheb;
	double* arrelsCheb;
	double* leg;
	double* arrelsLeg;
	double* coefsCheb;
	double* coefsLeg;
	double* nodesTrap;
	double tol = 1e-8;

	// Llegim el nombre de nodes amb els que integrarem
	printf("Introduïu el nombre de nodes: ");
	if(scanf("%d", &n) != 1) {
	  fprintf(stderr, "ERROR: Introduïu un nombre natural de nodes");
	  return 1;
	}
	
	// Inicialitzem les llistes
	cheb = (double *) malloc((n+1) * sizeof(double));
	arrelsCheb = (double *) malloc(n * sizeof(double));
	leg = (double *) malloc((n+1) * sizeof(double));
	arrelsLeg = (double *) malloc(n * sizeof(double));
	coefsCheb = (double *) malloc(n * sizeof(double));
	coefsLeg = (double *) malloc(n * sizeof(double));
	nodesTrap = (double *) malloc(n * sizeof(double));
	
	// Calculem les arrels del polinomi de Chebyshev i els coeficients per a la integral
	chebyshev(n, cheb);
	trobarIntervals(cheb, n, arrelsCheb);
	int i = 0;
	for(i = 0; i < n; i++) {
		arrelsCheb[i] = newton(cheb, n, arrelsCheb[i], tol);
	}
	coeficientsCheb(n, coefsCheb);

	// Calculem la integral
	double intCheb = 0;
	for(i = 0; i < n; i++) {
		intCheb += coefsCheb[i] * fabs(arrelsCheb[i]) * sqrt(1 - arrelsCheb[i] * arrelsCheb[i]);	
	}

	// Fem el mateix càlcul amb Gauss-Legendre
	legendre(n, leg);
	trobarIntervals(leg, n, arrelsLeg);
	for(i = 0; i < n; i++) {
		arrelsLeg[i] = newton(leg, n, arrelsLeg[i], tol);
	}
	coeficientsLeg(n, leg, arrelsLeg, coefsLeg);

	// Calculem la integral
	double intLeg = 0;
	for(i = 0; i < n; i++) {
		intLeg += coefsLeg[i] * fabs(arrelsLeg[i]);	
	}	

	// Calculem la integral amb la regla dels trapezis composta
	double d = 2./n;
	for(i = 0; i < n; i++) {
		nodesTrap[i] = fabs(-1 + i*d);	
	}
	double intTrap = trapezis(n, -1, 1, nodesTrap);

	printf("Quadratura per Gauss-Chebyshev = %lf\n", intCheb);
	printf("Quadratura per Gauss-Legendre = %lf\n", intLeg);
	printf("Quadratura per Trapezis = %lf\n", intTrap);
}
