#include "localdensityapproximation.h"
#include "system.h"
#include "Orbitals/contractedgaussian.h"
#include "Orbitals/gaussianprimitive.h"
#include <iostream>
#include <iomanip>

using arma::mat;

LocalDensityApproximation::LocalDensityApproximation(System*    system,
                                                     arma::mat* densityMatrix) :
        ExchangeCorrelationFunctional(system) {
    m_densityMatrix = densityMatrix;
}

double LocalDensityApproximation::evaluateEnergy(double rho) {
    const double rs  = pow(3.0/(4.0*3.1415926535897932384*rho), 1.0/3.0);
    return (rho < 1e-15) ? 0 : epsilonX(rs, rho) + epsilonC(rs, rho);
    //return epsilonX(rs,rho);
    //return epsilonC(rs,rho);
}

double LocalDensityApproximation::evaluatePotential(double rho) {
    const double rs  = pow(3.0/(4*3.1415926535897932384*rho), 1.0/3.0);
    return  (rho < 1e-20) ? 0 : vX(rs, rho) + vC(rs, rho);
    //return vC(rs,rho);
}

double LocalDensityApproximation::epsilonX(double rs, double rho) {
    double cx = -(3.0/4.0)*pow(3.0/3.1415926535897932384, 1.0/3.0);
    return cx * pow(rho, 1.0/3.0);
}

double LocalDensityApproximation::vX(double rs, double rho) {
    double cx = -(3.0/4.0)*pow(3.0/3.1415926535897932384, 1.0/3.0);
    return (4.0/3.0) * cx * pow(rho, 1.0/3.0);
}

double LocalDensityApproximation::vC(double rs, double rho) {
    const double x  = sqrt(rs);
    double       Ai  =  0.0621814 / 2.0;
    const double xi = -0.10498;
    const double bi =  3.72744;
    const double ci = 12.9352;
    const double Qi = sqrt(4*ci - bi*bi);
    const double Xi = xi*xi + bi*xi + ci;
    const double X  = x * x + bi*x  + ci;

    const double eps =  Ai * (log(x*x/X) + 2*bi/Qi*atan(Qi/(2*x+bi)) - bi*xi/Xi * (  log((x-xi)*(x-xi)/X) + 2*(bi+2*xi)/Qi*atan(Qi/(2*x+bi))));
    const double deps_drho = -Ai*rs/(6*x) * (   (2*ci+bi*x)/(x*(ci+x*(bi+x)))
                                             -(4*bi)/(Qi*Qi+(bi+2*x)*(bi+2*x))
                                             -(bi*xi)/(ci+xi*(bi+xi)) * (-(bi+2*x) / (ci+x*(bi+x))
                                                                         +2/(x-xi)
                                                                         -4*(bi+2*xi)/(Qi*Qi + (bi+2*x)*(bi+2*x))));
    return eps + deps_drho;
}

double LocalDensityApproximation::epsilonC(double rs, double rho) {
    ////
    //// http://theory.rutgers.edu/~giese/notes/DFT.pdf
    ////
    double x  = sqrt(rs);

    double A  =  0.0621814 / 2.0;
    double xi = -0.10498;
    double bi =  3.72744;
    double ci = 12.9352;
    double Qi = sqrt(4*ci - bi*bi);
    double Xi = xi*xi + bi*xi + ci;
    double X  = x * x + bi*x  + ci;

    return A * (   log(x*x/X) + 2*bi/Qi*atan(Qi/(2*x+bi)) - bi*xi/Xi * (  log((x-xi)*(x-xi)/X) + 2*(bi+2*xi)/Qi*atan(Qi/(2*x+bi))    )   );
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
    /*double x     = sqrt(rs);
    double rho43 = pow(rho, 4.0/3.0);
    return      (-0.982717+x*(-4.80211+x*(-4.12098+x*(4.47538+x*(0.51543+(0.00263132+0.00321452*x)*x)))))
            /   (x*(0.409286+x)*(0.409286+x)*(1824.98*x-85.4376*x*x*x + x*x*x*x*x)*rho43);
    */
    //3.1415926535897932384
    double x  = sqrt(rs);

    double A  =  0.0621814 / 2.0;
    double xi = -0.10498;
    double bi =  3.72744;
    double ci = 12.9352;
    double Qi = sqrt(4*ci - bi*bi);
    double Xi = xi*xi + bi*xi + ci;
    double X  = x * x + bi*x  + ci;

    double dcdx =  -A*(bi*x*xi*((Qi*Qi + (bi + 2*x)*(bi + 2*x))*(2*bi*x + 2*ci + 2*x*x -
                    (bi + 2*x)*(x - xi)) - 4*(bi + 2*xi)*(x - xi)*(bi*x + ci + x*x))
                    + 4*bi*x*(x - xi)*(bi*x + ci + x*x)*(bi*xi + ci + xi*xi) - (Qi*Qi
                    + (bi + 2*x)*(bi + 2*x))*(x - xi)*(bi*xi + ci + xi*xi)*(2*bi*x + 2*ci +
                    2*x*x - x*(bi + 2*x)))/(x*(Qi*Qi + (bi + 2*x)*(bi + 2*x))*(x - xi)*(bi*x
                    + ci + x*x)*(bi*xi + ci + xi*xi));
    double dxdrs = 1/(2*sqrt(rs));
    double drsdrho = -pow(6.0, 1.0/3.0)*pow(1.0/rho,1.0/3.0)/(6*pow(3.1415926535897932384,1.0/3.0)*rho);
    return dcdx * dxdrs * drsdrho;
}

double LocalDensityApproximation::dEpsilonX(double rs, double rho) {
    double rho23 = pow(rho, 2.0/3.0);
    return - 0.984745 / rho23;
    //double drsdrho  = 4*pow(2.0, 2.0/3.0)*pow(m_pi, 4.0/3.0) / (3*pow(3.0, 1.0/3.0) * pow(1.0/(rs*rs*rs), 4.0/3.0));
    //return drsdrho * (3*pow(3.0/(2*m_pi), 2.0/3.0) / (2*rs*rs));
}

