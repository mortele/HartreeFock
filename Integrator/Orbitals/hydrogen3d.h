#pragma once
#include "orbital.h"
#include <cmath>

class Hydrogen3D : public Orbital {
protected:
    double m_rMax       = 35.0;
    double m_thetaMax   = M_PI;
    double m_phiMax     = 2*M_PI;

    double solidHarmonics(int l, int m, double r, double theta, double phi);
    double normalization(int n, int l, int m);

public:
    Hydrogen3D() : Orbital(3,3) {}
    double  computeWavefunction(double *coordinates, int *quantumNumbers);
    double  integrandOne(double *allCoordinates, int *allQuantumNumbers);
    double  integrandTwo(double *allCoordinates, int *allQuantumNumbers);
    double* getCoordinateScales();
    void    updateCoordinateScales(int* allQuantumNumbers, int numberOfQuantumNumbers);
};

