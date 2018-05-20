#include<stdio.h>
#include<math.h>
#include"polinomis.h"

// Avalua el polinomi de Legendre de grau n a x
double legendre(double x, int n) {
  if(n < 0) {
    fprintf(stderr, "ERROR: El polinomi de Legendre només es pot calcular per n >= 0.");
		return -1;
  } else if(n == 0) {
    return 1;
  } else if(n == 1) {
    return x;
  }

  return ((2*n - 1)*x*legendre(x, n-1) - (n-1)*legendre(x, n-2)) / (double) n;
}

// Avalua el polinomi de Chebyshev de grau n a x
double chebyshev(double x, int n) {
  if(n < 0) {
    fprintf(stderr, "ERROR: El polinomi de Chebyshev només es pot calcular per n >= 0.");
		return -1;
  } else if(n == 0) {
    return 1;
  } else if(n == 1) {
    return x;
  }

  return 2*x*chebyshev(x, n-1) - chebyshev(x, n-2);
}

// Busca els punts on el polinomi de Legendre de grau n canvia de signe
void trobarIntervalsLeg(int n, double* arrels) {
	// Nombre de punts on avaluem
	int N = 20;
	double d;
	int nArrels = 0;

	int i;

	while(nArrels != n) {
		d = 2./N;
		nArrels = 0;
		for(i = 0; i < N; i++) {
			if(legendre(-1 + i*d, n) * legendre(-1 + (i+1)*d, n) < 0) {
				arrels[nArrels] = -1 + (2*i + 1)*d / 2.;
				nArrels++;
			}
		}
		N = 10*N;
	}
}

// Busca els punts on el polinomi de Chebyshev de grau n canvia de signe
void trobarIntervalsCheb(int n, double* arrels) {
	// Nombre de punts on avaluem
	int N = 20;
	double d;
	int nArrels = 0;

	int i;

	while(nArrels != n) {
		d = 2./N;
		nArrels = 0;
		for(i = 0; i < N; i++) {
			if(chebyshev(-1 + i*d, n) * chebyshev(-1 + (i+1)*d, n) < 0) {
				arrels[nArrels] = -1 + (2*i + 1)*d / 2.;
				nArrels++;
			}
		}
		N = 10*N;
	}
}

// Avalua la derivada del polinomi de Legendre de grau n a x
double legendreDerivada(double x, int n) {
	if(n < 0) {
		fprintf(stderr, "ERROR: El polinomi de Legendre només es pot calcular per n >= 0."); 
	} else if(n == 0) {
		return 0;
	}

	return (-n*x*legendre(x,n) + n*legendre(x,n-1)) / (1-(x*x));
}

// Avalua la derivada del polinomi de Chebyshev de grau n a x
double chebyshevDerivada(double x, int n) {
	if(n < 0) {
		fprintf(stderr, "ERROR: El polinomi de Chebyshev només es pot calcular per n >= 0."); 
	} else if(n == 0) {
		return 0;
	}

	return (-n*x*chebyshev(x,n) + n*chebyshev(x,n-1)) / (1-(x*x));
}

// Busca l'arrel del polinomi de Legendre de grau n propera a x amb tolerància tol
double newtonLeg(double x, int n, double tol) {
	// Hi anirem guardant els punts de la successió
	double xn = x;	

	double tempX;

	do {
		tempX = xn;
		xn = xn - legendre(xn,n)/legendreDerivada(xn,n);
	} while(fabs(xn - tempX) > tol);

	return xn;
}

// Busca l'arrel del polinomi de Chebyshev de grau n continguda a l'interval [a,b] amb tolerància tol
double newtonCheb(double x, int n, double tol) {
	// Hi anirem guardant els punts de la successió
	double xn = x;	

	double tempX;

	do {
		tempX = xn;
		xn = xn - chebyshev(xn,n)/chebyshevDerivada(xn,n);
	} while(fabs(xn - tempX) > tol);

	return xn;
}

