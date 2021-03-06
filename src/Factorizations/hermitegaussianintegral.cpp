#include "Factorizations/hermitegaussianintegral.h"

using arma::cube;
using arma::zeros;
using arma::vec;
using arma::zeros;
using std::max;
using std::cout;
using std::endl;


HermiteGaussianIntegral::HermiteGaussianIntegral() :
        m_t(0),
        m_u(0),
        m_v(0),
        m_nucleusPosition(zeros<vec>(3)),
        m_boysFunction(BoysFunction()){
}


double HermiteGaussianIntegral::getCoefficient(int n, int t, int u, int v) {
    if (t >= 0 && u >= 0 && v >= 0) {
        return m_coefficients(n)(t,u,v);
    } else {
        return 0;
    }
}


void HermiteGaussianIntegral::setupCoefficients(GaussianPrimitive* primitive1,
                                                GaussianPrimitive* primitive2,
                                                vec nucleusPosition) {
    int t = primitive1->xExponent() + primitive2->xExponent();
    int u = primitive1->yExponent() + primitive2->yExponent();
    int v = primitive1->zExponent() + primitive2->zExponent();

    m_nucleusPosition   = nucleusPosition;
    double  alpha       = primitive1->exponent();
    double  beta        = primitive2->exponent();
    double  p           = alpha + beta;
    vec     P           = (alpha * primitive1->nucleusPosition() +
                           beta  * primitive2->nucleusPosition()) / p;
    vec     PC          = P - nucleusPosition;
    setupCoefficients(t,u,v,p,PC);
}

void HermiteGaussianIntegral::setupCoefficients(int         t,
                                                int         u,
                                                int         v,
                                                double      p,
                                                arma::vec   PC) {
    m_PC         = PC;
    m_t          = t;
    m_u          = u;
    m_v          = v;
    m_tuv        = m_t + m_u + m_v;
    int maxIndex = m_tuv+1;
    m_coefficients.set_size(maxIndex);
    arma::field<arma::cube>& R = m_coefficients;
    BoysFunction& F = m_boysFunction;
    for (int i=0; i<maxIndex; i++) {
        m_coefficients(i) = zeros<cube>(maxIndex+1, maxIndex+1, maxIndex+1);
    }

    double x = p * arma::dot(m_PC, m_PC);
    m_coefficients(0)(0,0,0) = F.computeAndApplyDownwardRecurrence(x, m_tuv+1);

    //          n                        n          n
    // Set all R    values according to R    = (-2p)  F (p CP*CP)
    //          000                      000           n
    double minusTwoPPowerM = 1;
    for (int m = 1; m < m_tuv+1; m++) {
        minusTwoPPowerM *= (-2*p);
        m_coefficients(m)(0,0,0) = minusTwoPPowerM * F[m];
    }

    for (int tuv = 1; tuv < m_tuv+1;     tuv++)
    for (int n   = 0; n   < m_tuv+1-tuv; n  ++)
    for (int t   = 0; t   < m_t+1;       t  ++)
    for (int u   = 0; u   < m_u+1;       u  ++)
    for (int v   = 0; v   < m_v+1;       v  ++) {
        if (t + u + v != tuv  ||  t + u + v == 0) {
            continue;
        }
        int tuvMax = max(t, max(u, v));

        double newCoefficient = 0;
        if (tuvMax == t) {
            newCoefficient = (t-1)    * getCoefficient(n+1,t-2,u,v) +
                             m_PC(0)  * getCoefficient(n+1,t-1,u,v);
        } else if (tuvMax == u) {
            newCoefficient = (u-1)    * getCoefficient(n+1,t,u-2,v) +
                             m_PC(1)  * getCoefficient(n+1,t,u-1,v);
        } else if (tuvMax == v) {
            newCoefficient = (v-1)    * getCoefficient(n+1,t,u,v-2) +
                             m_PC(2)  * getCoefficient(n+1,t,u,v-1);
        }
        m_coefficients(n)(t,u,v) = newCoefficient;
    }
}


























