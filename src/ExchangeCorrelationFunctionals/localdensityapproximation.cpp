#include "localdensityapproximation.h"
#include "system.h"
#include "Orbitals/contractedgaussian.h"
#include <iostream>

using arma::mat;

LocalDensityApproximation::LocalDensityApproximation(System*    system,
                                                     arma::mat* densityMatrix) :
        ExchangeCorrelationFunctional(system) {
    m_densityMatrix = densityMatrix;
    m_A   =  0.0621814;
    m_x0  = -0.409286;
    m_b   = 13.0720;
    m_c   = 42.7198;
    m_pi  = 3.1415926535897932384;
    m_Q   = sqrt(4*m_c - m_b*m_b);
    m_Xx0 = m_x0*m_x0 + m_b*m_x0 + m_c;
}

double LocalDensityApproximation::evaluate(double x, double y, double z, int p, int q) {
    m_A   =  0.0621814;
    m_x0  = -0.409286;
    m_b   = 13.0720;
    m_c   = 42.7198;
    m_pi  = 3.1415926535897932384;
    m_Q   = sqrt(4*m_c - m_b*m_b);
    m_Xx0 = m_x0*m_x0 + m_b*m_x0 + m_c;

    const arma::mat& P = *m_densityMatrix;
    ContractedGaussian* Gp = m_system->getBasis().at(p);
    ContractedGaussian* Gq = m_system->getBasis().at(q);
    const double rho = P(p,q) * Gp->evaluate(x,y,z) * Gq->evaluate(x,y,z);
    const double rs  = pow(3.0/(4*m_pi*rho), 1.0/3.0);
    return (rho==0 ? 0 : epsilonx(rs) + epsilonc(rs));
}

double LocalDensityApproximation::epsilonc(double rs) {
    m_A   =  0.0621814;
    m_x0  = -0.409286;
    m_b   = 13.0720;
    m_c   = 42.7198;
    m_pi  = 3.1415926535897932384;
    m_Q   = sqrt(4*m_c - m_b*m_b);
    m_Xx0 = m_x0*m_x0 + m_b*m_x0 + m_c;

    double x  = sqrt(rs);
    double Xx = X(x);
    return  m_A/2.0 * ( log(x*x/Xx) + (2*m_b/m_Q)*atan(m_Q/(2*x+m_b)) -m_b*m_x0/m_Xx0
                        *(log((x-m_x0)*(x-m_x0)/Xx) + 2*(m_b + 2*m_x0)/m_Q * atan(m_Q/(2*x+m_b))));
}

double LocalDensityApproximation::epsilonx(double rs) {
    return -3.0*pow(9.0/(32*m_pi*m_pi),1.0/3.0) * (1.0/rs);
}
