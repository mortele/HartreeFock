#pragma once
#include "Orbitals/orbital.h"
#include <cmath>

class Hydrogen3D : public Orbital {
protected:
    double m_rMax       = 5.0;
    double m_thetaMax   = M_PI;
    double m_phiMax     = 2*M_PI;

public:
    Hydrogen3D() : Orbital(3,3) {}
    double  computeWavefunction(double *coordinates, int *numberOfQuantumNumbers);
    double* getCoordinateScales();
};

