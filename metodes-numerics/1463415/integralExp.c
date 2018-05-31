#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"polinomis.h"
#include"matrius.h"

double sqexp(double);

int main() {
	// Fem el càlcul de la integral de exp(-x^2)/sqrt(1-x^2) entre -1 i 1 amb Gauss-Chebyshev
	int n;
	double* cheb;
	double* arrelsCheb;
	double* leg;
	double* arrelsLeg;
	double** ACheb;
	double* bCheb;
	double* coefsCheb;
	double** ALeg;
	double* bLeg;
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
	ALeg = (double**) malloc(n * sizeof(double*));
	ACheb = (double**) malloc(n * sizeof(double*));
	bLeg = (double*) malloc(n * sizeof(double*));
	bCheb = (double*) malloc(n * sizeof(double*));

	int i,j;
	for(i = 0; i < n; i++) {
		ALeg[i] = (double*) malloc((n+1) * sizeof(double));	
		ACheb[i] = (double*) malloc((n+1) * sizeof(double));	
	}

	// Calculem les arrels del polinomi de Chebyshev
	chebyshev(n, cheb);
	trobarIntervals(cheb, n, arrelsCheb);
	for(i = 0; i < n; i++) {
		arrelsCheb[i] = newton(cheb, n, arrelsCheb[i], tol);
	}

	// Fem el mateix càlcul amb Gauss-Legendre
	legendre(n, leg);
	trobarIntervals(leg, n, arrelsLeg);
	for(i = 0; i < n; i++) {
		arrelsLeg[i] = newton(leg, n, arrelsLeg[i], tol);
	}

	// Omplim les matrius de Vandermonde i les resolem
	for(i = 0; i < n; i++) {
		ALeg[0][i] = 1;	
		ACheb[0][i] = 1;	
		for(j = 1; j < n; j++) {
			ALeg[j][i] = ALeg[j - 1][i] * arrelsLeg[i];
			ACheb[j][i] = ACheb[j - 1][i] * arrelsCheb[i];
		}
	}

	bLeg[0] = 2;
	bCheb[0] = M_PI;
	for(i = 1; i < n; i++) {
		if(i % 2 == 1) {
			bCheb[i] = 0;
			bLeg[i] = 0;
		} else {
			bLeg[i] = 2./(i+1);
			bCheb[i] = (double)(i-1)/i * bCheb[i-2];
		}
	}

	triangula(n, ALeg, bLeg);
	triangula(n, ACheb, bCheb);

	subsEndarrera(n, ALeg, bLeg, coefsLeg);
	subsEndarrera(n, ACheb, bCheb, coefsCheb);


	// Calculem la integral per Chebyshev
	double intCheb = 0;
	for(i = 0; i < n; i++) {
		intCheb += coefsCheb[i] * exp(-arrelsCheb[i] * arrelsCheb[i]);	
	}

	// Calculem la integral per Legendre
	double intLeg = 0;
	for(i = 0; i < n; i++) {
	  intLeg += coefsLeg[i] * sqexp(arrelsLeg[i]); 
	}

	// Calculem la integral amb la regla dels trapezis composta
	double d = 2./n;
	for(i = 0; i < n; i++) {
		nodesTrap[i] = sqexp(-1 + i*d);	
	}
	double intTrap = trapezis(n, -1 + 0.1, 1 - 0.1, nodesTrap);

	printf("Quadratura per Gauss-Chebyshev = %lf\n", intCheb);
	printf("Quadratura per Gauss-Legendre = %lf\n", intLeg);
}

double sqexp(double x) {
	return exp(-x*x)/sqrt(1 - x*x);
}
