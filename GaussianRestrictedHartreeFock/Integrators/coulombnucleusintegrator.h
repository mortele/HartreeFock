#pragma once
#include <armadillo>
#include <iostream>
#include "hermitegaussian.h"


class CoulombNucleusIntegrator {
private:
    arma::vec       m_nucleusPosition;
    HermiteGaussian m_hermiteGaussian;

public:
    CoulombNucleusIntegrator();
    void setNucleusPosition(arma::vec nucleusPosition);
};

