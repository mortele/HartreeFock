#pragma once
#include "Orbitals/gaussianprimitive.h"
#include "Factorizations/hermitegaussian.h"

class OverlapIntegrator {
private:
    double          m_Ex                = 0;
    double          m_Ey                = 0;
    double          m_Ez                = 0;
    double          m_sqrtPiOverP       = 0;
    HermiteGaussian m_hermiteGaussian;


public:
    OverlapIntegrator();
    double computeIntegral(GaussianPrimitive* primitive1,
                           GaussianPrimitive* primitive2);
    double getIntegralIndicesDimension(int i, int j, int dimension);
    double getIntegralX() { return m_Ex; }
    double getIntegralY() { return m_Ey; }
    double getIntegralZ() { return m_Ez; }
    double getIntegralDimension(int dimension);
};

