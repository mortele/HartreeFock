//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <cassert>
#include <boost/timer.hpp>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Solvers/unrestrictedhartreefock.h"
#include "Atoms/atom.h"
#include "Atoms/hydrogen.h"
#include "Atoms/oxygen.h"
#include "Atoms/helium.h"
#include "examples.h"

#include "Integrators/overlapintegrator.h"
#include "Integrators/kineticintegrator.h"
#include "Integrators/electronnucleusintegrator.h"
#include "Orbitals/gaussianprimitive.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {
    //Examples::Hm();
    //Examples::He();
    //Examples::HeHp();
    //Examples::H2();
    //Examples::H20();

    vec nucleus1  {0,            0,             0};
    vec nucleus2  {0,            0,             1};

    System system = System(2);
    system.addAtom(new Oxygen  ("6-311++G**", vec{0,0,0}));
    system.addAtom(new Hydrogen("6-311++G**", vec{1.809,0,0}));

    cout << system.twoBodyElements(10,21,8,26) << endl;


    return 0;
}
