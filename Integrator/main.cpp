#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include "Orbitals/orbital.h"
#include "Orbitals/harmonicoscillator2d.h"
#include "montecarlointegrator.h"
#include "ran1.h"

using namespace std;
long         idum      = -time(0);
const int    n         = (int) 1e7;
const double rmax      = 3.0;
const double eps       = pow(10,-12);
const double factorOne = (2*M_PI) * rmax;
const double factorTwo = factorOne*factorOne;
const double reference = 1.253314137; // 1111
//const double reference = 0.939985603; // 1212
//const double reference = 0.3133285343; // 1221
//const double reference = 0.744155269; // 1616
//const double reference = 0.3035370176; // 2346


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

double psi(double r, double theta, int n, int m) {
    double x = r*r;
    double oneOverSquareRoot = 1/sqrt(M_PI*factorial(n+fabs(m)));
    double rPowerM = (m==0) ? 1 : ((m==1)||(m==-1) ? r : r);
    return L(x, n, m) * oneOverSquareRoot * exp(-x/2.0) * rPowerM;
}

double psi1(double r, double theta) {
    return psi(r, theta, 0, 0);
}
double psi2(double r, double theta) {
    return psi(r, theta, 0, -1);
}
double psi3(double r, double theta) {
    return psi(r, theta, 0, +1);
}
double psi4(double r, double theta) {
    return psi(r, theta, 0, -2);
}

double integrandOne(double r1,
                    double theta1) {
    return psi1(r1, theta1) * psi1(r1, theta1) * r1;
}

double integrandTwo(double r1,
                    double theta1,
                    double r2,
                    double theta2) {
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    int m4 = 0;
    double r12 = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(theta2-theta1));
    double phase = cos((m3-m1)*theta1)*cos((m4-m2)*theta2) - sin((m3-m1)*theta1)*sin((m4-m2)*theta2);
    return (r12 < eps) ? 0 : phase * psi1(r1, theta1) * psi1(r2, theta2) * psi1(r1, theta1) * psi1(r2, theta2) * r1 * r2 / r12;
}

void hei(int* input) {
    cout << input[0] << ", " << input[1] << endl;
}

int* mapToOrbitals(int p) {
    int* quantumNumbers = new int[2];
    switch (p) {
        case 1:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 0;
            break;
        case 2:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = -1;
            break;
        case 3:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 1;
            break;
        case 4:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = -2;
            break;
        case 5:
            quantumNumbers[0] = 1;
            quantumNumbers[1] = 0;
            break;
        case 6:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 2;
            break;
        default:
            cout << "Invalid orbital <" << p << ">." << endl;
            break;
    }
    return quantumNumbers;
}

int* generateQuantumNumbersOne(int* indices) {
    int* allQuantumNumbers = new int[4];
    for (int i=0; i<2; i++) {
        int* quantumNumbers = mapToOrbitals(indices[i]);
        allQuantumNumbers[2*i+0] = quantumNumbers[0];
        allQuantumNumbers[2*i+1] = quantumNumbers[1];
    }
    return allQuantumNumbers;
}

int* generateQuantumNumbersTwo(int* indices) {
    int* allQuantumNumbers = new int[8];
    for (int i=0; i<4; i++) {
        int* quantumNumbers = mapToOrbitals(indices[i]);
        allQuantumNumbers[2*i+0] = quantumNumbers[0];
        allQuantumNumbers[2*i+1] = quantumNumbers[1];
    }
    return allQuantumNumbers;
}

int main() {

    //int indices [] = {1,1};
    //int* allQuantumNumbers = generateAllQuantumNumbersOne(indices);

    int indices [] = {1,1,1,1};
    int* allQuantumNumbers = generateQuantumNumbersTwo(indices);


    MonteCarloIntegrator* MCInt = new MonteCarloIntegrator();
    MCInt->setOrbital(new HarmonicOscillator2D());
    //cout << "Integral: " << MCInt->integrateOne(allQuantumNumbers) << endl;
    cout << "Integral: " << MCInt->integrateTwo(allQuantumNumbers) << endl;
    cout << "stdDev:   " << MCInt->getStandardDeviation() << endl;


    // Timing of the code.
    clock_t start, finish;
    start = clock();

    // Integration loop.
    double sum  = 0;
    double sum2 = 0;
    double term = 0;
    for(int i = 0; i < n; i++) {
        double r1     = ran1(&idum)*rmax;
        double r2     = ran1(&idum)*rmax;
        double theta1 = ran1(&idum)*2*M_PI;
        double theta2 = ran1(&idum)*2*M_PI;
        term  = integrandTwo(r1, theta1, r2, theta2) * factorTwo;
        //term  = integrandOne(r1, theta1) * factorOne;
        sum  += term;
        sum2 += term*term;
    }
    sum  /= ((double) n);
    sum2 /= ((double) n*n);

    // Terminal output.
    finish = clock();
    cout << endl << " * Importance sampled Monte Carlo " << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "Value of numerical approx. I' = " << (sum) << endl;
    cout << "Exact value of integral I     = " << reference << endl;
    cout << "Time usage                    = " << (finish-start)/1000000.0;
    cout << " [seconds]" << endl;
    cout << "Variance                      = " << ((sum2) - (sum)*(sum))/n << endl;
    cout << "Standard deviation            = " << sqrt(((sum2) - (sum)*(sum))/n) << endl;
    cout << "Relative error                = " << fabs(reference-(sum))/reference << endl;
    cout << "Number of points, n           = " << "10^" << log10(n) << endl << endl;
    return 0;
}
