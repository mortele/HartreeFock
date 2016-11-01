#include "unrestrictedhartreefock.h"

using std::cout;
using std::endl;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;


UnrestrictedHartreeFock::UnrestrictedHartreeFock(System* system) :
        HartreeFock(system) {

    m_numberOfSpinUpElectrons = 0;
    m_numberOfSpinDownElectrons = 0;

    m_epsilonUp                   = zeros(m_numberOfBasisFunctions);
    m_epsilonDown                 = zeros(m_numberOfBasisFunctions);
    m_epsilonOldUp                = zeros(m_numberOfBasisFunctions);
    m_epsilonOldDown              = zeros(m_numberOfBasisFunctions);
    m_fockMatrixUp                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixDown              = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeUp           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeDown         = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixUp         = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixDown       = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixTildeUp    = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixTildeDown  = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_densityMatrixUp             = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_densityMatrixDown           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);

}
