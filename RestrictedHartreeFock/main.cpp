#include <iostream>
#include <restrictedhartreefock.h>

using namespace std;

int main()
{

    int nrOfParticles = 2;
    int nrOfSpinOrbitals = 12;
    RestrictedHartreeFock* rhf = new RestrictedHartreeFock(nrOfParticles,nrOfSpinOrbitals);
    rhf->computeSolutionBySCF();
    return 0;
}

