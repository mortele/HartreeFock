#pragma once
#include <armadillo>
#include <iostream>
#include "Factorizations/hermitegaussian.h"
#include "Factorizations/hermitegaussianintegral.h"


class CoulombNucleusIntegrator {
private:
    arma::vec               m_nucleusPosition;
    HermiteGaussian         m_hermiteGaussian;
    HermiteGaussianIntegral m_hermiteGaussianIntegral;

public:
    CoulombNucleusIntegrator();
    void setNucleusPosition(arma::vec nucleusPosition);
    double computeIntegral(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2);
    double computeIntegral(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2,
                           arma::vec          nucleusPosition);
};

