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

public:
    HermiteGaussianIntegral();
    void setupCoefficients(GaussianPrimitive* primitive1,
                           GaussianPrimitive* primitive2,
                           arma::vec          nucleusPosition);
    void setupCoefficients(int t, int u, int v, double p,
                           arma::vec PC);
    double getCoefficient(int n, int t, int u, int v);
};

