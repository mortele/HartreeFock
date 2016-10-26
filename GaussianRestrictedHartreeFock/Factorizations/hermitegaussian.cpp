#include "Factorizations/hermitegaussian.h"
#include <cmath>
#include <iostream>
#include <iomanip>


using std::max;
using std::exp;
using std::cout;
using std::endl;
using arma::zeros;
using arma::cube;
using arma::vec;

HermiteGaussian::HermiteGaussian() {
}

void HermiteGaussian::setupCoefficients(GaussianPrimitive& primitive1,
                          GaussianPrimitive& primitive2) {

    m_nucleusPosition1  = primitive1.nucleusPosition();
    m_nucleusPosition2  = primitive2.nucleusPosition();
    m_xExponent1        = primitive1.xExponent();
    m_yExponent1        = primitive1.yExponent();
    m_zExponent1        = primitive1.zExponent();
    m_xExponent2        = primitive2.xExponent();
    m_yExponent2        = primitive2.yExponent();
    m_zExponent2        = primitive2.zExponent();
    m_xMaximumAngularMomentum = max(m_xExponent1, m_xExponent2);
    m_yMaximumAngularMomentum = max(m_yExponent1, m_yExponent2);
    m_zMaximumAngularMomentum = max(m_zExponent1, m_zExponent2);
    m_maximumAngularMomentum  = max(max(max(m_xMaximumAngularMomentum,
                                            m_yMaximumAngularMomentum),
                                        m_zMaximumAngularMomentum),
                                    1);
    m_exponent1         = primitive1.exponent();
    m_exponent2         = primitive2.exponent();
    m_exponentSum       = m_exponent1 + m_exponent2;

    m_coefficients[0] = zeros<cube>(m_xMaximumAngularMomentum+1,
                                    m_xMaximumAngularMomentum+1,
                                    2 * m_xMaximumAngularMomentum + 2);
    m_coefficients[1] = zeros<cube>(m_yMaximumAngularMomentum+1,
                                    m_yMaximumAngularMomentum+1,
                                    2 * m_yMaximumAngularMomentum + 2);
    m_coefficients[2] = zeros<cube>(m_zMaximumAngularMomentum+1,
                                    m_zMaximumAngularMomentum+1,
                                    2 * m_zMaximumAngularMomentum + 2);
    computeCoefficients();
}

double HermiteGaussian::getCoefficientDimension(int i, int j, int dimension) {
    if (dimension == 0) {
        return getCoefficientX(i,j);
    } else if (dimension == 1) {
        return getCoefficientY(i,j);
    } else if (dimension == 2) {
        return getCoefficientZ(i,j);
    }
}

double HermiteGaussian::getCoefficientDimension(int i, int j, int k, int dimension) {
    return m_coefficients[dimension](i,j,k);
}

bool HermiteGaussian::isCoefficientNonZero(int i, int j, int t) {
    if (t < 0 || t > (i+j) || i < 0 || j < 0) {
        return false;
    } else {
        return true;
    }
}


void HermiteGaussian::computeCoefficients() {
    int iA_loopLimits [] = {m_xExponent1 + 1, m_yExponent1 + 1, m_zExponent1 + 1};
    int iB_loopLimits [] = {m_xExponent2 + 1, m_yExponent2 + 1, m_zExponent2 + 1};
    int t_loopLimits  [] = {m_xExponent1 + m_xExponent2 + 1,
                            m_yExponent1 + m_yExponent2 + 1,
                            m_zExponent1 + m_zExponent2 + 1};

    double  alpha   = m_exponent1;
    double  beta    = m_exponent2;
    double  p       = alpha + beta;
    double  mu      = alpha * beta / p;
    vec     AB      = m_nucleusPosition1 - m_nucleusPosition2;
    vec     P       = (alpha * m_nucleusPosition1 + beta * m_nucleusPosition2) / p;
    vec     PA      = P - m_nucleusPosition1;
    vec     PB      = P - m_nucleusPosition2;

    for (int i = 0; i<3; i++) {
        cube& E = m_coefficients[i];
        double AB_ = AB(i);
        double PA_ = PA(i);
        double PB_ = PB(i);

        int iA = 0;
        // Need to place E(0,0,0) here to ensure it is set for loopLimits = (0,0,1).
        E(0,0,0) = exp(- mu * AB_ * AB_);
        for (int iB = 0; iB < iB_loopLimits[i]; iB++) {
            for (int t = 0; t < t_loopLimits[i]; t++) {
                if (iA == 0 && iB == 0 && t == 0) {
                    //E(0,0,0) = exp(- mu * AB_ * AB_);
                } else {
                    // E(i, j-1, t-1)
                    double previousIBpreviousT = 0;
                    if (isCoefficientNonZero(iA, iB-1, t-1)) {
                        previousIBpreviousT = E(iA, iB-1, t-1);
                    }

                    // E(i, j-1, t)
                    double previousIB = 0;
                    if (isCoefficientNonZero(iA, iB-1, t)) {
                        previousIB = E(iA, iB-1, t);
                    }

                    // E(i,j-1, t+1)
                    double previousIBnextT = 0;
                    if (isCoefficientNonZero(iA, iB-1, t+1)) {
                        previousIBnextT = E(iA, iB-1, t+1);
                    }
                    // E(i,j,t) = E(i,j-1,t-1)/2p + Xpb E(i,j-1,t) + (t+1) E(i,j-1,t+1)
                    E(iA,iB,t) = (1./(2*p))   * previousIBpreviousT   +
                                 PB_          * previousIB            +
                                 (t+1)        * previousIBnextT;
                }
            }
        }

        for (iA = 1; iA < iA_loopLimits[i]; iA++) {
            for (int iB = 0; iB < iB_loopLimits[i]; iB++) {
                for (int t = 0; t < t_loopLimits[i]; t++) {
                    // E(i-1, j, t-1)
                    double previousIApreviousT = 0;
                    if (isCoefficientNonZero(iA-1, iB, t-1)) {
                        previousIApreviousT = E(iA-1, iB, t-1);
                    }

                    // E(i-1, j, t)
                    double previousIA = 0;
                    if (isCoefficientNonZero(iA-1, iB, t)) {
                        previousIA = E(iA-1, iB, t);
                    }

                    // E(i-1,j, t+1)
                    double previousIAnextT = 0;
                    if (isCoefficientNonZero(iA-1, iB, t+1)) {
                        previousIAnextT = E(iA-1, iB, t+1);
                    }

                    // E(i,j,t) = E(i-1,j,t-1)/2p + Xpb E(i-1,j,t) + (t+1) E(i-1,j,t+1)
                    E(iA,iB,t) = (1./(2*p))   * previousIApreviousT   +
                                 PA_          * previousIA            +
                                 (t+1)        * previousIAnextT;
                }
            }
        }
    }
}

