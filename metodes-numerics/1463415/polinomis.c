#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"polinomis.h"

// Avalua el polinomi de grau 'n' amb coeficients emmagatzemats a 'p' al punt 'x'
double avalua(double* p, int n, double x) {
  int i;
  double resultat = p[n];

  for(i = n - 1; i >= 0; i--) {
    resultat = p[i] + resultat*x;
  }

	return resultat;
}

// Guarda a 'der' els coeficients de la derivada del polinomi de grau 'n' amb coeficients emmagatzemats a 'p'
void deriva(double* p, int n, double* der) {
  int i;
  for(i = 0; i < n; i++) {
    der[i] = (i+1)*p[i+1];
  }
}

// Guarda els coeficients de p+q a la llista r
void suma(double* p, int n, double* q, int m, double* r) {
  int i;
  if(n == m) {
    for(i = 0; i <= n; i++) {
      r[i] = p[i] + q[i];
    }
  } else if(n < m) {
    for(i = 0; i <= n; i++) {
      r[i] = p[i] + q[i];
    }
    for(i = n + 1; i <= m; i++) {
      r[i] = q[i];
    }
  } else if(m < n) {
    for(i = 0; i <= m; i++) {
      r[i] = p[i] + q[i];
    }
    for(i = m + 1; i <= n; i++) {
      r[i] = p[i];
    }
  }
}

// Multplica el polinomi 'p' per l'escalar 'a' i ho guarda a 'ap'
void escala(double a, double* p, int n, double* ap) {
  int i;
  for(i = 0; i <= n; i++) {
    ap[i] = a*p[i];
  }
}

// Guarda a la llista 'L' els coeficients del polinomi de Legendre de grau 'n'
void legendre(int n, double* L) {
  if(n == 0) {
    L[0] = 1;
  } else if (n == 1) {
    L[0] = 0;
    L[1] = 1;
  } else {
		int i,j;
 		double** tempL = (double**) malloc((n+1) * sizeof(double*)); 
		for(i = 0; i <= n; i++) {
			tempL[i] = (double*) malloc((n+1) * sizeof(double));
		}
		double* tempL1 = (double*) malloc((n+1) * sizeof(double));
		double* tempL2 = (double*) malloc((n+1) * sizeof(double));


		// Inicialitzem la iteració: L[0] serà el polinomi de Legendre de grau 0, i L[1] el de grau 1
		tempL[0][0] = 1;
		tempL[0][n] = 0;

		tempL[1][0] = 0;
		tempL[1][1] = 1;

		// Realitzem la iteració
		for(i = 2; i <= n ; i++) {
			// Guardem els dos anteriors polinomis de Legendre
			for(j = 0; j <= n; j++) {
				tempL1[j] = tempL[i-1][j];
				tempL2[j] = tempL[i-2][j];
			}

			// Multipliquem tempL1 per (2i - 1)/i
			escala((double)(2*i - 1)/i, tempL1, i-1, tempL1);
			// Multipliquem tempL1 per x	
			for(j = i; j > 0 ; j--) {
				tempL1[j] = tempL1[j - 1];	
			}
			tempL1[0] = 0;
			// Multipliquem tempL2 per (1 - i)/i
			escala((double)(1-i)/i, tempL2, i-2, tempL2);
			// Sumem i ho guardem a tempL
			suma(tempL1, i, tempL2, i-2, tempL[i]);
		}

		// Guardem el resultat a L
		for(i = 0; i <= n; i++) {
			L[i] = tempL[n][i];
		}
	}
}

// Guarda a la llista 'C' els coeficients del polinomi de Chebyshev de grau 'n'
void chebyshev(int n, double* C){
  if(n == 0) {
    C[0] = 1;
  } else if (n == 1) {
    C[0] = 0;
    C[1] = 1;
  }

	int i,j;
 	double** tempC = (double**) malloc((n+1) * sizeof(double*)); 
	for(i = 0; i <= n; i++) {
		tempC[i] = (double*) malloc((n+1) * sizeof(double));
	}
	double* tempC1 = (double*) malloc((n+1) * sizeof(double));
	double* tempC2 = (double*) malloc((n+1) * sizeof(double));


	// Inicialitzem la iteració: C[0] serà el polinomi de Chebyshev de grau 0, i C[1] el de grau 1
	tempC[0][0] = 1;
	tempC[0][n] = 0;

	tempC[1][0] = 0;
	tempC[1][1] = 1;

	// Realitzem la iteració
	for(i = 2; i <= n ; i++) {
		// Guardem els dos anteriors polinomis de Chebyshev
		for(j = 0; j <= n; j++) {
			tempC1[j] = tempC[i-1][j];
			tempC2[j] = tempC[i-2][j];
		}
		
		// Multipliquem tempC1 per 2
		escala(2, tempC1, i-1, tempC1);
		// Multipliquem tempC1 per x	
		for(j = i; j > 0 ; j--) {
			tempC1[j] = tempC1[j - 1];	
		}
		tempC1[0] = 0;
		// Multipliquem tempC2 per -1
		escala(-1, tempC2, i-2, tempC2);
		// Sumem i ho guardem a tempC
		suma(tempC1, i, tempC2, i-2, tempC[i]);
	}

	// Guardem el resultat a C
	for(i = 0; i <= n; i++) {
		C[i] = tempC[n][i];
	}
}

// Busca els punts on el polinomi de Legendre o Chebyshev de grau n canvia de signe
void trobarIntervals(double* p, int n, double* arrels) {
	// Nombre de punts on avaluem
	int N = 11;
	double d;
	int nArrels = 0;

	int i;

	while(nArrels != n) {
		d = 2./N;
		nArrels = 0;
		for(i = 0; i < N; i++) {
			if((avalua(p, n, -1 + i*d) * avalua(p, n, -1 + (i+1)*d) < 0)) {
				arrels[nArrels] = -1 + (2*i + 1)*d / 2.;
				nArrels++;
			}
		}
		N = 10*N;
	}
}

// Busca l'arrel del polinomi de Legendre o Chebyshev de grau n propera a x amb tolerància tol
double newton(double* p, int n, double x, double tol) {
	// Hi anirem guardant els punts de la successió
	double xn = x;	

	double tempX;

	double* derP = (double*) malloc((n+1) * sizeof(double));
	deriva(p, n, derP);

	do {
		tempX = xn;
		xn = xn - avalua(p, n, xn)/avalua(derP, n, xn);
	} while(fabs(xn - tempX) > tol);

	return xn;
}

// Omple la llista 'coefs' amb els coeficients de la quadratura de Gauss-Legendre de grau n
void coeficientsLeg(int n, double* L, double* arrels, double* coefs) {
	double* derL = (double*) malloc((n+1) * sizeof(double));
	deriva(L, n, derL);

	int i;
	for(i = 0; i < n; i++) {
		coefs[i] = 2 / ( (1-arrels[i]*arrels[i]) * (avalua(derL, n, arrels[i])*avalua(derL, n, arrels[i])) );
	}
}

// Omple la llista 'coefs' amb els coeficients de la qudratura de Gauss-Chebyshev de grau n
void coeficientsCheb(int n, double* coefs) {
	int i;
	for(i = 0; i < n; i++) {
		coefs[i] = M_PI / (double) n;
	}
}

// Calcula la integral de exp(-x^2)/sqrt(1 - x^2) a l'interval [a,b] partint-lo en n - 1 intervals
double trapezis(int n, double a, double b, double* f) {
	double h = (b - a)/(n - 1);

	double s = f[0] + f[n - 1];
	int i;
	for(i = 1; i <= n - 2; i++) {
		s += 2*f[i];
	}

	return (h/2.) * s;
}

