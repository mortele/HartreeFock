#pragma once
#include <armadillo>


class GaussianPrimitive {
    friend class ContractedGaussian;
    friend class GaussianPrimitiveProduct;

private:
    int         m_xExponent                 = 0;
    int         m_yExponent                 = 0;
    int         m_zExponent                 = 0;
    int         m_maximumAngularMomentum    = 0;
    double      m_exponent                  = 0;
    double      m_constantTerm              = 1.0;
    arma::vec   m_nucleusPosition           = arma::zeros<arma::vec>(3);
    double      m_coefficient               = 0;


public:
    GaussianPrimitive(int i, int j, int k, double a, arma::vec nucleusPosition, double constantTerm=1.0);
    void setCoefficient(double coefficient);
    double evaluate(arma::vec& r);
    static GaussianPrimitive product(GaussianPrimitive& primitive1, GaussianPrimitive& primitive2);

    double      exponent()                  const { return m_exponent;  }
    int         xExponent()                 const { return m_xExponent; }
    int         yExponent()                 const { return m_yExponent; }
    int         zExponent()                 const { return m_zExponent; }
    int         maximumAngularMomentum()    const { return m_maximumAngularMomentum; }
    arma::vec   nucleusPosition()           const { return m_nucleusPosition; }

};

