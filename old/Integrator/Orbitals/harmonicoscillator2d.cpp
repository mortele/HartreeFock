#include <cmath>
#include <iostream>
#include "harmonicoscillator2d.h"

using std::cout;
using std::endl;

double HarmonicOscillator2D::computeWavefunction(double* coordinates,
                                                 int*    quantumNumbers) {
    double  r       = coordinates[0];
    int     n       = quantumNumbers[0];
    int     m       = quantumNumbers[1];
    int     mAbs    = m < 0 ? -m : m;

    double x                 = r * r;
    double L                 = Orbital::associatedLaguerrePolynomial(x, n, mAbs);
    double rPowerM           = std::pow(r, mAbs);
    double exponential       = std::exp(-0.5*x);
    double oneOverSquareRoot = 1. / std::sqrt(M_PI * Orbital::factorial(n + mAbs));
    return L * rPowerM * exponential * oneOverSquareRoot;
}

double HarmonicOscillator2D::integrandOne(double* allCoordinates,
                                          int*    allQuantumNumbers) {

    double r                    = allCoordinates[0];
    double theta                = allCoordinates[1];
    int    m1                   = allQuantumNumbers[1];
    int    m2                   = allQuantumNumbers[3];
    double integrationMeasure   = r;

    double phase         = std::cos(theta * (m1 - m2));
    double waveFunction1 = computeWavefunction(allCoordinates, allQuantumNumbers);
    double waveFunction2 = computeWavefunction(allCoordinates, allQuantumNumbers+2);
    double integrand     = integrationMeasure * waveFunction1 * waveFunction2 * phase;
    return integrand;
}

double HarmonicOscillator2D::integrandTwo(double* allCoordinates,
                                          int*    allQuantumNumbers) {
    double r1                   = allCoordinates[0];
    double theta1               = allCoordinates[1];
    double r2                   = allCoordinates[2];
    double theta2               = allCoordinates[3];
    int    m1                   = allQuantumNumbers[1];
    int    m2                   = allQuantumNumbers[3];
    int    m3                   = allQuantumNumbers[5];
    int    m4                   = allQuantumNumbers[7];
    double integrationMeasure   = r1 * r2;

    double r12           = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(theta2-theta1));
    double oneOverR12    = r12 < 1e-12 ? 0 : 1./r12;
    double phase         = cos((m3-m1)*theta1)*cos((m4-m2)*theta2) -
                           sin((m3-m1)*theta1)*sin((m4-m2)*theta2);
    double waveFunction1 = computeWavefunction(allCoordinates,   allQuantumNumbers  );
    double waveFunction2 = computeWavefunction(allCoordinates+2, allQuantumNumbers+2);
    double waveFunction3 = computeWavefunction(allCoordinates,   allQuantumNumbers+4);
    double waveFunction4 = computeWavefunction(allCoordinates+2, allQuantumNumbers+6);
    double integrand     = integrationMeasure * oneOverR12 * phase *
                           waveFunction1 * waveFunction2 * waveFunction3 * waveFunction4;
    return integrand;
}

double* HarmonicOscillator2D::getCoordinateScales() {
    double* scales = new double[m_dimensions];
    scales[0] = m_rMax;
    scales[1] = m_thetaMax;
    return scales;
}
