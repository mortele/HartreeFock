#include "examples.h"
#include <armadillo>
#include <restrictedhartreefock.h>

using std::cout;
using std::endl;

int Examples::twoDimensionalQuantumDot() {
    int nrOfParticles = 2;
    int nrOfSpinOrbitals = 12;

    arma::vec oneBodyElements = arma::zeros<arma::vec>(nrOfSpinOrbitals/2);

    oneBodyElements(0) = 1.0; oneBodyElements(1) = 2.0; oneBodyElements(2) = 2.0;
    oneBodyElements(3) = 3.0; oneBodyElements(4) = 3.0; oneBodyElements(5) = 3.0;

    RestrictedHartreeFock* rhf = new RestrictedHartreeFock(nrOfParticles,nrOfSpinOrbitals);
    rhf->setAnalyticOneBodyElements(oneBodyElements);
    bool tableLoaded = rhf->setIntegralTable("../Integrator/IntegralTables/HO_2d_6_nonzero.dat");

    //cout << tableLoaded << endl;

    if (tableLoaded) {
        rhf->computeSolutionBySCF();
    }

    return 0;

}
