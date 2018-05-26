#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"polinomis.h"

int main() {
	int n = 8;
	double arrels[8];
	double tol = 0.000001;

	double* C = (double*) malloc((n + 1) * sizeof(double));
	chebyshev(n, C);
	trobarIntervals(C, n, arrels);

	int i = 0;
	double d = 2/11.;

	for(i = 0; i < n; i++) {
		arrels[i] = newton(C, n, arrels[i], tol);
		printf("%lf\n", arrels[i]);
	}

	return 0;
}
