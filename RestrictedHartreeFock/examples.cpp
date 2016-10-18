#include "examples.h"
#include <armadillo>
#include <restrictedhartreefock.h>

using std::cout;
using std::endl;

int Examples::twoDimensionalQuantumDot() {

    int nrOfParticles = 20;
    int nrOfSpinOrbitals = 20;

    arma::vec oneBodyElements = arma::zeros<arma::vec>(nrOfSpinOrbitals/2);

    for (int i=0; i<nrOfSpinOrbitals/2; i++) {
        if (i==0) {
            oneBodyElements(i) = 1.0;
        } else if (i==1 || i==2) {
            oneBodyElements(i) = 2.0;
        } else if (i==3 || i==4 || i==5) {
            oneBodyElements(i) = 3.0;
        } else if (i==6 || i==7 || i==8 || i==9) {
            oneBodyElements(i) = 4.0;
        }
    }

    /*oneBodyElements(0) = 1.0; oneBodyElements(1) = 2.0; oneBodyElements(2) = 2.0;
    oneBodyElements(3) = 3.0; oneBodyElements(4) = 3.0; oneBodyElements(5) = 3.0;
    oneBodyElements(6) = 4.0; oneBodyElements(7) = 4.0; oneBodyElements(8) = 4.0;
    oneBodyElements(9) = 4.0; */

    RestrictedHartreeFock* rhf = new RestrictedHartreeFock(nrOfParticles,nrOfSpinOrbitals);
    rhf->setAnalyticOneBodyElements(oneBodyElements);
    bool tableLoaded = rhf->setIntegralTable("../Integrator/IntegralTables/test2.dat");

    //cout << tableLoaded << endl;

    if (tableLoaded) {
        rhf->computeSolutionBySCF();
    }

    return 0;

}
