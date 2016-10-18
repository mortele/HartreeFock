#pragma once
#include "Orbitals/orbital.h"
#include "Orbitals/hydrogen3d.h"
#include <cmath>

class HarmonicOscillator2D : public Orbital {
protected:
    double m_rMax     = 5.0;
    double m_thetaMax = 2*M_PI;

public:
    HarmonicOscillator2D() : Orbital(2,2) {}
    double  computeWavefunction (double* coordinates, int* quantumNumbers);
    double  integrandOne(double* allCoordinates, int* allQuantumNumbers);
    double  integrandTwo(double* allCoordinates, int* allQuantumNumbers);
    double* getCoordinateScales();
};
