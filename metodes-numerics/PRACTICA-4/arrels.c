#include<stdio.h>
#include<math.h>
#include"polinomis.h"

int main() {
	int n = 4;
	double arrels[4];
	double tol = 0.000001;

	trobarIntervalsLeg(n, arrels);

	int i = 0;
	for(i = 0; i < n; i++) {
		arrels[i] = newtonLeg(arrels[i], n, tol);
		printf("arrel %d = %lf\n", i, arrels[i]);
	}

	return 0;
}
