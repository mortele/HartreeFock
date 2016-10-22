#pragma once
#include <armadillo>
#include "gaussianprimitive.h"

class HermiteGaussian {
private:
    int         m_xExponent1, m_xExponent2  = 0;
    int         m_yExponent1, m_yExponent2  = 0;
    int         m_zExponent1, m_zExponent2  = 0;
    double      m_exponent1,  m_exponent2   = 0;
    double      m_exponentSum               = 0;
    int         m_xMaximumAngularMomentum   = 0;
    int         m_yMaximumAngularMomentum   = 0;
    int         m_zMaximumAngularMomentum   = 0;
    int         m_maximumAngularMomentum    = 0;

    arma::vec   m_nucleusPosition1;
    arma::vec   m_nucleusPosition2;
    arma::cube  m_coefficients[3];

    void    computeCoefficients();
    bool    isCoefficientNonZero(int t, int i, int j);

public:
    HermiteGaussian();
    void    set(GaussianPrimitive& primitive1, GaussianPrimitive& primitive2);
    double  getCoefficientX(int i, int j) const { return m_coefficients[0](i,j,0); }
    double  getCoefficientY(int i, int j) const { return m_coefficients[1](i,j,0); }
    double  getCoefficientZ(int i, int j) const { return m_coefficients[2](i,j,0); }
    double  getCoefficientDimension(int i, int j, int dimension);
    double  getExponentSum()              const { return m_exponentSum; }
};
