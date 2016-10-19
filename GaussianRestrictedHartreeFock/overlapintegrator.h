#pragma once
#include "gaussianprimitive.h"
#include "hermitegaussian.h"

class OverlapIntegrator {
private:
    HermiteGaussian m_hermiteGaussian;

public:
    OverlapIntegrator();
    double computeIntegral(GaussianPrimitive& primitive1,
                           GaussianPrimitive& primitive2);

};

