#pragma once
#include "Orbitals/orbital.h"
#include <cmath>

class HarmonicOscillator2D : public Orbital {
private:
    double m_rMax     = 3.0;
    double m_thetaMax = 2*M_PI;

    double associatedLaguerrePolynomial(double x, int n, int m);
    int factorial(int n);

public:
    HarmonicOscillator2D() : Orbital(2,2) {}
    double computeWavefunction (double* coordinates, int* quantumNumbers);
    double integrandOne(double* allCoordinates, int* allQuantumNumbers);
    double integrandTwo(double *allCoordinates, int *allQuantumNumbers);
    double* getCoordinateScales();
};
