#include<stdio.h>
#include<math.h>

// Avalua el polinomi de Legendre de grau n a x
double legendre(double x, int n) {
    if(n < 0) {
        fprintf(stderr, "ERROR: El polinomi de Legendre només es pot calcular per n >= 0.");
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
        fprintf(stderr, "ERROR: El polinomi de Legendre només es pot calcular per n >= 0.");
    } else if(n == 0) {
        return 1;
    } else if(n == 1) {
        return x;
    }
    
    return 2*x*chebyshev(x, n-1) - chebyshev(x, n-2);
}
