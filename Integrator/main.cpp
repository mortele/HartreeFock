#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include "ran1.h"

using namespace std;
long         idum      = -time(0);
const int    n         = (int) 5e7;
const double rmax      = 3.0;
const double eps       = pow(10,-12);
const double factor    = (2*M_PI) * (2*M_PI) * rmax * rmax;
const double reference = 1.253314137;


double psi1(double r) {
    return exp(-r*r/2.0)/sqrt(M_PI);
}
double psi2(double r) {
    return exp(-r*r/2.0)/sqrt(M_PI);
}
double psi3(double r) {
   return exp(-r*r/2.0)/sqrt(M_PI);
}
double psi4(double r) {
    return exp(-r*r/2.0)/sqrt(M_PI);
}

double integrand(double r1,
                 double theta1,
                 double r2,
                 double theta2) {

    double r12 = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(theta2-theta1));
    return (r12 < eps) ? 0 : psi1(r1) * psi2(r2) * psi3(r1) * psi4(r2) * r1 * r2 / r12;
}

int main() {
    // Timing of the code.
    clock_t start, finish;
    start = clock();

    // Integration loop.
    double term, sum = 0, sum2 = 0;
    for(int i = 0; i < n; i++) {
        float r1     = ran1(&idum)*rmax;
        float r2     = ran1(&idum)*rmax;
        float theta1 = ran1(&idum)*2*M_PI;
        float theta2 = ran1(&idum)*2*M_PI;
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
    cout << "Value of numerical approx. I' = " << sum << endl;
    cout << "Exact value of integral I     = " << reference << endl;
    cout << "Time usage                    = " << (finish-start)/1000000.0;
    cout << " [seconds]" << endl;
    cout << "Variance                      = " << (sum2 - sum*sum)/n << endl;
    cout << "Standard deviation            = " << sqrt((sum2 - sum*sum)/n) << endl;
    cout << "Relative error                = " << fabs(reference-sum)/reference << endl;
    cout << "Number of points, n           = " << "10^" << log10(n) << endl << endl;
    return 0;
}
