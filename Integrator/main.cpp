#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <string>
#include "Orbitals/orbital.h"
#include "Orbitals/harmonicoscillator2d.h"
#include "Orbitals/hydrogen3d.h"
#include "montecarlointegrator.h"
#include "ran1.h"
#include "integraltable.h"

using std::cout;
using std::endl;

int main() {
    IntegralTable table;
    MonteCarloIntegrator integrator;
    integrator.setOrbital(new HarmonicOscillator2D());
    std::string fileName = "../Integrator/IntegralTables/test2.dat";
    table.createTwoBodyTable(fileName, 0, 10, (int) 1e7, integrator);

    return 0;

    for (int i=0; i<5; i++) {
        int qm [] = {i,i};
        int* allQm = Orbital::generateQuantumNumbers(qm, 1, 1);
        double I = integrator.integrateOne(allQm, (int) 5e7);
        cout << "i:" << i << ", I=" << I << endl;
    }

    return 0;
}
