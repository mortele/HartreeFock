#pragma once
#include <armadillo>
#include "Orbitals/gaussianprimitive.h"
#include "Math/boysfunction.h"


class HermiteGaussianIntegral {
private:
    int                         m_t;
    int                         m_u;
    int                         m_v;
    int                         m_tuv;
    int                         m_maxExponents;
    arma::vec                   m_PC;
    arma::vec                   m_nucleusPosition;
    arma::field<arma::cube>     m_coefficients;
    BoysFunction                m_boysFunction;

    int computeMaximumExponents(GaussianPrimitive& primitive1,
                                GaussianPrimitive& primitive2);
    double getCoefficient(int n, int t, int u, int v);

public:
    HermiteGaussianIntegral();
    void setupCoefficients(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2,
                           arma::vec          nucleusPosition);
};

