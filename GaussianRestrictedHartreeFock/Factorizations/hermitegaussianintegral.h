#pragma once
#include <armadillo>
#include "Orbitals/gaussianprimitive.h"

class HermiteGaussianIntegral {
private:
    int         m_t;
    int         m_u;
    int         m_v;
    arma::vec   m_nucleusPosition;
    arma::cube  m_coefficients;

public:
    HermiteGaussianIntegral();
    void setupCoefficients(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2,
                           arma::vec          nucleusPosition);
};

