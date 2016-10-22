#include "Integrators/coulombnucleusintegrator.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

CoulombNucleusIntegrator::CoulombNucleusIntegrator() :
        m_nucleusPosition(zeros<vec>(3)),
        m_hermiteGaussian(HermiteGaussian()) {

}

void CoulombNucleusIntegrator::setNucleusPosition(arma::vec nucleusPosition) {
    m_nucleusPosition = nucleusPosition;
}
