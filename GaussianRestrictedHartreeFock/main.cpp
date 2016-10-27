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
#include "Atoms/hydrogen.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::pow;
using arma::vec;
using arma::zeros;

int main(int, char**) {

    vec nucleus1 = {0, 0, 0};
    //vec nucleus2 = {0, 0, 1.4};

    Hydrogen H1(nucleus1);
    //Hydrogen H2(nucleus2);




    return 0;
}

