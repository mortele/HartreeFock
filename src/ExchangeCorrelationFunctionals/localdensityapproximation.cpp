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

double LocalDensityApproximation::evaluateEnergy(double x, double y, double z, int p, int q) {
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
    return (rho<1e-20 ? 0 : epsilonX(rs,rho) + epsilonC(rs,rho));
}

double LocalDensityApproximation::evaluatePotential(double x, double y, double z, int p, int q) {

    double pi  = 3.1415926535897932384;


    const arma::mat& P = *m_densityMatrix;
    ContractedGaussian* Gp = m_system->getBasis().at(p);
    ContractedGaussian* Gq = m_system->getBasis().at(q);
    const double rho = P(p,q) * Gp->evaluate(x,y,z) * Gq->evaluate(x,y,z);
    const double rs  = pow(3.0/(4*pi*rho), 1.0/3.0);
    //return (rho<1e-20 ? 0 : epsilonX(rs,rho) + epsilonC(rs,rho) + rho*(dEpsilonC(rs,rho) + dEpsilonX(rs,rho)));
    return (rho<1e-20 ? 0 : epsilonX(rs,rho) + epsilonC(rs,rho));
}

double LocalDensityApproximation::epsilonX(double rs, double rho) {
    //double pi  = 3.1415926535897932384;
    //return -3.0*pow(9.0/(32*pi*pi),1.0/3.0) * (1.0/rs);
    return - 0.4582 / rs;
}

double LocalDensityApproximation::epsilonC(double rs, double rho) {
    //m_A   =  0.0621814;
    //m_x0  = -0.409286;
    //m_b   = 13.0720;
    //m_c   = 42.7198;
    //m_pi  = 3.1415926535897932384;
    //m_Q   = sqrt(4*m_c - m_b*m_b);
    //m_Xx0 = m_x0*m_x0 + m_b*m_x0 + m_c;

    //double x  = sqrt(rs);
    //double Xx = X(x);
    //return  m_A/2.0 * ( log(x/Xx) + (2*m_b/m_Q)*atan(m_Q/(2*x+m_b)) -m_b*m_x0/m_Xx0
    //                    *(log((x-m_x0)*(x-m_x0)/Xx) + 2*(m_b + 2*m_x0)/m_Q * atan(m_Q/(2*x+m_b))));
    double x = sqrt(rs);
    return      18.10330000 * atan(0.0224499 / (-6.536 + x))
            +    2.41869000 * atan(0.0224499 / ( 6.536 + x))
            +    0.03109070 * log (x / (42.7198+x*(13.072+x)))
            +    0.00443137 * log ((0.409286+x)*(0.409286+x) / (42.7198+x*(13.072+x)));
}

double LocalDensityApproximation::dEpsilonC(double rs, double rho) {
    //m_A   =  0.0621814;
    //m_x0  = -0.409286;
    //m_b   = 13.0720;
    //m_c   = 42.7198;
    //m_pi  = 3.1415926535897932384;
    //m_Q   = sqrt(4*m_c - m_b*m_b);
    //m_Xx0 = m_x0*m_x0 + m_b*m_x0 + m_c;
    //
    //double x = sqrt(rs);
    //double aa = 0.14809777061418503; // (3/4pi)^(4/3)
    //double rho43 = /* aa */ pow(rho, 4.0/3.0);
    //return (-0.0460076+(-0.238897-0.309044*x)*x)/(x*x*(0.409286+x)*(0.409286+x)*(42.7198+x*(13.072+x))*rho43);
    //double mathematica = (-0.0230038+x*(-0.119448+x*(-0.136785+(0.0446515+0.00321452*x)*x))) /
    //                      (x*x*(0.409286+x)*(0.409286+x)*(42.7198+x*(13.072+x))*rho43);
    //
    // //double x  = sqrt(rs);
    //double Xx = X(x);
    //
    //double drsdrho  = -(1.0/rho43)/(pow(6.0, 2.0/3.0) * pow(m_pi,1.0/3.0)); //4*pow(2.0, 2.0/3.0)*pow(m_pi, 4.0/3.0) / (3*pow(3.0, 1.0/3.0) * pow(1.0/(rs*rs*rs), 4.0/3.0));
    //double dxdrs    = 1.0/(2*x);
    //
    // // ln (x² / X)
    //double firstTerm    = (m_c-x*x) / (x*x*(m_b+x) + m_c*x);
    //
    // // (2b / Q) tan⁻¹ (Q / (2x + b))
    //double secondTerm   = - 4*m_b / ((m_b+2*x)*(m_b+2*x) + m_Q*m_Q);
    //
    // // ln ((x-x0)² / X(x))
    //double thirdTerm    = (m_b*x + m_x0*(m_b+2*x) + 2*m_c) / ((x-m_x0) * Xx);
    //
    // // (2(b+2x0) / Q) tan⁻¹ (Q / (2x + b))
    //double fourthTerm   = - 4*(m_b+2*m_x0) / (m_b*m_b + 4*m_b*x + m_Q*m_Q + 4*x*x);
    //
    //double mine = drsdrho * dxdrs * m_A/2.0 * (firstTerm + secondTerm + (-m_b*m_x0/m_Xx0) * (thirdTerm + fourthTerm));
    // //return drsdrho * dxdrs * m_A/2.0 * (firstTerm + secondTerm + (-m_b*m_x0/m_Xx0) * (thirdTerm + fourthTerm));
    //
    // //std::cout << std::fabs(firstTerm-mathFirst) << std::endl;
    //return mathematica;
    double x     = sqrt(rs);
    double rho43 = pow(rho, 4.0/3.0);
    return      (-0.982717+x*(-4.80211+x*(-4.12098+x*(4.47538+x*(0.51543+(0.00263132+0.00321452*x)*x)))))
            /   (x*(0.409286+x)*(0.409286+x)*(1824.98*x-85.4376*x*x*x + x*x*x*x*x)*rho43);

}

double LocalDensityApproximation::dEpsilonX(double rs, double rho) {
    double rho23 = pow(rho, 2.0/3.0);
    return - 0.984745 / rho23;
    //double drsdrho  = 4*pow(2.0, 2.0/3.0)*pow(m_pi, 4.0/3.0) / (3*pow(3.0, 1.0/3.0) * pow(1.0/(rs*rs*rs), 4.0/3.0));
    //return drsdrho * (3*pow(3.0/(2*m_pi), 2.0/3.0) / (2*rs*rs));
}

