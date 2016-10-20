#pragma once
#include <armadillo>
#include "overlapintegrator.h"
#include "gaussianprimitive.h"


class KineticIntegrator {
private:
    arma::vec           m_overlapIntegrals;
    arma::vec           m_T;
    arma::mat           m_adjustedOverlapIntegrals;
    OverlapIntegrator   m_overlapIntegrator;
    GaussianPrimitive   m_primitive1;
    GaussianPrimitive   m_primitive2;

    void computeAdjustedOverlapIntegral(int dimension, int adjustment);
    void computeT(int dimension);

public:
    KineticIntegrator();
    double computeIntegral(GaussianPrimitive& primitive1, GaussianPrimitive& primitive2);

};
