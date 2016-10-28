#include "restrictedhartreefock.h"

using std::cout;
using std::endl;
using arma::zeros;
using arma::mat;
using arma::vec;


RestrictedHartreeFock::RestrictedHartreeFock(System* system) :
        m_system(system) {
    m_numberOfBasisFunctions    = m_system->getNumberOfBasisFunctions();
    m_numberOfElectrons         = m_system->getAtoms().size();
    m_basis                     = m_system->getBasis();
    m_epsilon                   = zeros<vec>(m_numberOfBasisFunctions);
    m_epsilonOld                = zeros<vec>(m_numberOfBasisFunctions);
    m_fockMatrix                = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_U                         = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_densityMatrix             = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
}


void RestrictedHartreeFock::computeFockMatrix() {

}


double RestrictedHartreeFock::solve(double convergenceCriteria) {

}
