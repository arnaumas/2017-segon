#include<stdio.h>
#include<math.h>
#include"polinomis.h"

int main() {
	// Fem el càlcul de la integral de exp(-x^2)/sqrt(1-x^2) entre -1 i 1 amb Gauss-Chebyshev
	int n = 8;
	double arrelsCheb[8];
	double coefsCheb[8];
	double tol = 1e-6;

	// Calculem les arrels del polinomi de Chebyshev i els coeficients per a la integral
	trobarIntervalsCheb(n, arrelsCheb);
	int i = 0;
	for(i = 0; i < n; i++) {
		arrelsCheb[i] = newtonCheb(arrelsCheb[i], n, tol);
	}
	coeficientsCheb(n, coefsCheb);

	// Calculem la integral
	double intCheb = 0;
	for(i = 0; i < n; i++) {
		intCheb += coefsCheb[i] * exp(-arrelsCheb[i] * arrelsCheb[i]);	
	}

	printf("%lf\n", intCheb);

}
