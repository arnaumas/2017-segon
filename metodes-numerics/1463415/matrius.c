#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Resol el sistema Ax = b mitjançant substitució endarrera fent servir que A és triangular superior, i guarda el resultat a x
void subsEndarrera(int n, double** A, double* b, double* x) {
	int i, j;

	double s;
	for(i = n - 1; i >= 0; i--) {
		s = 0;
		for(j = i + 1; j < n; j++) {
			s += A[i][j]*x[j];	
		}
		x[i] = (b[i] - s)/A[i][i];
	}

}

// Triangula inferiorment la matriu associada al sistema Ax = b
void triangula(int n, double** A, double*b) {
	int i, j, k;
	double a = 0;

	for(i = 0; i < n; i++) {
		for(j = i + 1; j < n; j++) {
			a = A[j][i]/A[i][i];
			for(k = i; k < n; k++) {
				A[j][k] = A[j][k] - a*A[i][k];	
			}
			b[j] = b[j] - a*b[i];
		}
	}
}
	
