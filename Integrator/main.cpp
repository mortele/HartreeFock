#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <complex>
#include "ran1.h"

using namespace std;
long         idum      = -time(0);
const int    n         = (int) 5e7;
const double rmax      = 3.0;
const double eps       = pow(10,-12);
const double factor    = (2*M_PI) * (2*M_PI) * rmax * rmax;
//const double reference = 1.253314137; // 1111
//const double reference = 0.939985603; // 1212
//const double reference = 0.3133285343; // 1221
//const double reference = 0.744155269; // 1616
const double reference = 0.3035370176; // 2346


double L(double x, int n, int m) {
    return n == 0 ? 1 : 1 - x + fabs(m);
}

double factorial(int n) {
    double fac = 1;
    for (int i=1; i<n+1; i++) {
        fac *= i;
    }
    return fac;
}

complex<double> psi(double r, double theta, int n, int m) {
    complex<double> i(0,1);
    double x = r*r;
    double oneOverSquareRoot = 1/sqrt(M_PI*factorial(n+fabs(m)));
    double rPowerM = (m==0) ? 1 : ((m==1)||(m==-1) ? r : x);
    complex<double> phase = exp(i*((double)m*theta));
    return L(x, n, m) * oneOverSquareRoot * exp(-x/2.0) * rPowerM * phase;
}

// (n=0, m=0)
complex<double> psi1(double r, double theta) {
    return psi(r, theta, 0, 0);
}

// (n=0, m=-1)
complex<double> psi2(double r, double theta) {
    return psi(r, theta, 0, -1);
}

// (n=0, m=+1)
complex<double> psi3(double r, double theta) {
   return psi(r, theta, 0, +1);
}

// (n=0, m=-2)
complex<double> psi4(double r, double theta) {
    return psi(r, theta, 0, -2);
}

// (n=1, m=0)
complex<double> psi5(double r, double theta) {
    return psi(r, theta, 1, 0);
}

// (n=0, m=+2)
complex<double> psi6(double r, double theta) {
    return psi(r, theta, 0, +2);
}

complex<double> integrand(double r1,
                 double theta1,
                 double r2,
                 double theta2) {

    double r12 = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(theta2-theta1));
    return (r12 < eps) ? 0 : conj(psi2(r1, theta1) * psi3(r2, theta2)) * psi4(r1, theta1) * psi6(r2, theta2) * r1 * r2 / r12;
}

int main() {
    // Timing of the code.
    clock_t start, finish;
    start = clock();

    // Integration loop.
    complex<double> sum  = 0;
    complex<double> sum2 = 0;
    complex<double> term = 0;
    for(int i = 0; i < n; i++) {
        double r1     = ran1(&idum)*rmax;
        double r2     = ran1(&idum)*rmax;
        double theta1 = ran1(&idum)*2*M_PI;
        double theta2 = ran1(&idum)*2*M_PI;
        term  = integrand(r1, theta1, r2, theta2);
        sum  += term;
        sum2 += term*term;
    }
    sum  *= factor          / n;
    sum2 *= factor * factor / n;

    // Terminal output.
    finish = clock();
    cout << endl << " * Importance sampled Monte Carlo " << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "Value of numerical approx. I' = " << real(sum) << endl;
    cout << "Imaginary contamination of I' = " << imag(sum) << endl;
    cout << "Exact value of integral I     = " << reference << endl;
    cout << "Time usage                    = " << (finish-start)/1000000.0;
    cout << " [seconds]" << endl;
    cout << "Variance                      = " << (real(sum2) - real(sum)*real(sum))/n << endl;
    cout << "Standard deviation            = " << sqrt((real(sum2) - real(sum)*real(sum))/n) << endl;
    cout << "Relative error                = " << fabs(reference-real(sum))/reference << endl;
    cout << "Number of points, n           = " << "10^" << log10(n) << endl << endl;
    return 0;
}
