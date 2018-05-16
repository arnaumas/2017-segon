#include<stdio.h>
#include<math.h>

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
void trobarIntervalsLeg(int n, double* intervals) {
	d = 0.1;
	int nArrels = 0;

	int i;

	while(nArrels + 1 != n) {
		nArrels = 0;
		intervals[0] = -1;
		for(i = 0; i < n; i++) {
			if(legendre(-1 + i*d, n) * legendre(intervals[nArrels], n) < 0) {
				nArrels++;
				intervals[nArrels] = -1 + i*d;
			}
		}
		d = d/10.;
	}
}

// Busca els punts on el polinomi de Chebyshev de grau n canvia de signe
void trobarIntervalsCheb(int n, double* intervals) {
	d = 0.1;
	int nArrels = 0;

	int i;

	while(nArrels + 1 != n) {
		nArrels = 0;
		intervals[0] = -1;
		for(i = 0; i < n; i++) {
			if(chebyshev(-1 + i*d, n) * chebyshev(intervals[nArrels], n) < 0) {
				nArrels++;
				intervals[nArrels] = -1 + i*d;
			}
		}
		d = d/10.;
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

// Busca l'arrel del polinomi de Legendre de grau n continguda a l'interval [a,b] amb tolerància tol
double newtonLeg(double a, double b, int n, double tol) {
	// Hi anirem guardant els punts de la successió
	double x = (b - a)/2;	
	double tempX;

	do {
		tempX = x;
		x = x - legendre(x,n)/legendreDerivada(x,n);
	} while(fabs(x - tempX) < tol);

	return x;
}

// Busca l'arrel del polinomi de Chebyshev de grau n continguda a l'interval [a,b] amb tolerància tol
double newtonCheb(double a, double b, int n, double tol) {
	// Hi anirem guardant els punts de la successió
	double x = (b - a)/2;	
	double tempX;

	do {
		tempX = x;
		x = x - chebyshev(x,n)/chebyshevDerivada(x,n);
	} while(fabs(x - tempX) < tol);

	return x;
}

