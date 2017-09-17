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

    IntegralTable           table;
    MonteCarloIntegrator    integrator;
    std::string             fileName    = "../IntegralTables/Hydrogen_3d_6.dat";
    integrator.setOrbital(new Hydrogen3D());
    table.createTwoBodyTable(fileName, 1, 6, (int) 1e7, integrator);

    /*
    for (int i=9; i<10; i++) {
        int qm [] = {i,i};
        int* allQm = Orbital::generateQuantumNumbers(qm, 1, 1);
        integrator.updateCoordinateScales(allQm, 6);
        double I = integrator.integrateOne(allQm, (int) 1e8);
        cout << "i:" << i << ", I=" << std::setprecision(10) <<  I
             << ",  std.dev: " << integrator.getStandardDeviation()
             << "     norm:    " << 1/std::sqrt(I) << endl;
    }*/

    return 0;
}
