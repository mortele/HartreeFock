#include "localdensityapproximation.h"
#include "system.h"
#include "Orbitals/contractedgaussian.h"

using arma::mat;

LocalDensityApproximation::LocalDensityApproximation(System*    system,
                                                     arma::mat* densityMatrix) :
        ExchangeCorrelationFunctional(system) {
    m_densityMatrix = densityMatrix;
}



double LocalDensityApproximation::evaluate(double x, double y, double z, int p, int q) {
    const arma::mat& P = *m_densityMatrix;
    ContractedGaussian* Gp = m_system->getBasis().at(p);
    ContractedGaussian* Gq = m_system->getBasis().at(q);
    const double rho = P(p,q) * Gp->evaluate(x,y,z) * Gq->evaluate(x,y,z);
    return -Cx() * pow(rho, 1.0/3.0);
}
