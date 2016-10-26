#pragma once
#include "Orbitals/gaussianprimitive.h"
#include "Factorizations/hermitegaussian.h"
#include "Factorizations/hermitegaussianintegral.h"
#include <armadillo>


class ElectronElectronIntegrator {
private:
    double                  m_2sqrtPiToThe5;
    HermiteGaussian         m_hermiteGaussian12;
    HermiteGaussian         m_hermiteGaussian34;
    HermiteGaussianIntegral m_hermiteGaussianIntegral;

    void setupHermiteGaussianIntegral(GaussianPrimitive& primitive1,
                                      GaussianPrimitive& primitive2,
                                      GaussianPrimitive& primitive3,
                                      GaussianPrimitive& primitive4);

public:
    ElectronElectronIntegrator();
    double computeIntegral(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2,
                           GaussianPrimitive& primitive3,
                           GaussianPrimitive& primitive4);

};

