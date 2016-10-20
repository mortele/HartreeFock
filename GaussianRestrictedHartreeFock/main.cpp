#include <iostream>
#include "gaussianprimitive.h"
#include "overlapintegrator.h"
#include "kineticintegrator.h"
#include "hermitegaussian.h"

using std::cout;
using std::endl;
using arma::vec;
using arma::zeros;

int main(int, char**) {

    vec A = {1,0,0};
    vec B = {0,0,1};

    int x1 = 0;
    int y1 = 1;
    int z1 = 0;

    int x2 = 2;
    int y2 = 1;
    int z2 = 0;

    GaussianPrimitive primitive1 = GaussianPrimitive(x1,y1,z1, 2.0, A);
    GaussianPrimitive primitive2 = GaussianPrimitive(x2,y2,z2, 13.0, B);

    //OverlapIntegrator integrator;
    KineticIntegrator integrator;
    cout << integrator.computeIntegral(primitive1, primitive2) << endl;



    return 0;
}

