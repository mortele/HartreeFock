#include <iostream>
#include <iomanip>
#include <time.h>
#include "Orbitals/gaussianprimitive.h"
#include "Integrators/overlapintegrator.h"
#include "Integrators/kineticintegrator.h"
#include "Integrators/electronnucleusintegrator.h"
#include "Integrators/electronelectronintegrator.h"
#include "Factorizations/hermitegaussian.h"
#include "Math/boysfunction.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::pow;
using arma::vec;
using arma::zeros;


int main(int, char**) {

    vec A       = {0,0,0};
    vec B       = {0,0,1};
    vec C       = {1,0,1};
    vec D       = {0,0,1};

    int x1 = 2;
    int y1 = 1;
    int z1 = 6;

    int x2 = 2;
    int y2 = 3;
    int z2 = 3;

    int x3 = 2;
    int y3 = 2;
    int z3 = 0;

    int x4 = 1;
    int y4 = 0;
    int z4 = 3;

    GaussianPrimitive primitive1 = GaussianPrimitive(x1,y1,z1, 2.0, A);
    GaussianPrimitive primitive2 = GaussianPrimitive(x2,y2,z2, 1.0, B);
    GaussianPrimitive primitive3 = GaussianPrimitive(x3,y3,z3, 3.0, C);
    GaussianPrimitive primitive4 = GaussianPrimitive(x4,y4,z4, 0.1, D);

    //OverlapIntegrator integrator;
    //KineticIntegrator integrator;
    //ElectronNucleusIntegrator integrator;
    ElectronElectronIntegrator integrator;
    cout << setprecision(15) << integrator.computeIntegral(primitive1,
                                                           primitive2,
                                                           primitive3,
                                                           primitive4) << endl;




    return 0;
}

