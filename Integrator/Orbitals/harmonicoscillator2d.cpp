#include <cmath>
#include <iostream>
#include "harmonicoscillator2d.h"

using std::cout;
using std::endl;

double HarmonicOscillator2D::computeWavefunction(double* coordinates,
                                                 int*    quantumNumbers) {
    double  r       = coordinates[0];
    //double  theta   = coordinates[1];
    int     n       = quantumNumbers[0];
    int     m       = quantumNumbers[1];
    int     mAbs    = m < 0 ? -m : m;

    double x                 = r * r;
    double L                 = associatedLaguerrePolynomial(x, n, m);
    double rPowerM           = std::pow(r, mAbs);
    double exponential       = std::exp(-0.5*x);
    double oneOverSquareRoot = 1. / std::sqrt(M_PI * factorial(n + mAbs));
    return L * rPowerM * exponential * oneOverSquareRoot;
}

double HarmonicOscillator2D::integrandOne(double* allCoordinates,
                                          int*    allQuantumNumbers) {

    double r                    = allCoordinates[0];
    double theta                = allCoordinates[1];
    //int    n1                   = allQuantumNumbers[0];
    int    m1                   = allQuantumNumbers[1];
    //int    n2                   = allQuantumNumbers[2];
    int    m2                   = allQuantumNumbers[3];
    double integrationMeasure   = r;

    double phase        = std::cos(theta * (m1-m2));
    double waveFunction = computeWavefunction(allCoordinates, allQuantumNumbers);
    double integrand    =  integrationMeasure * waveFunction * waveFunction * phase;
    return integrand;
}

double* HarmonicOscillator2D::getCoordinateScales() {
    double* scales = new double[m_dimensions];
    scales[0] = m_rMax;
    scales[1] = m_thetaMax;
    return scales;
}

double HarmonicOscillator2D::associatedLaguerrePolynomial(double x,
                                                          int    n,
                                                          int    m) {
    return (n==0) ? 1 : (1 - x + std::fabs(m));
}

int HarmonicOscillator2D::factorial(int n) {
    int factorial = 1;
    for (int i=1; i<n+1; i++) {
        factorial *= i;
    }
    return factorial;
}
