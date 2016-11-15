//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <cassert>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Solvers/unrestrictedhartreefock.h"
#include "Atoms/atom.h"
#include "Atoms/hydrogen.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {

    vec nucleus1  {0, 0, 0};
    vec nucleus2  {0, 0, 1};
    vec nucleus3  {0, 1, 0};
    vec nucleus4  {1, 0, 0};
    vec nucleus5  {1, 1, 0};
    vec nucleus6  {0, 1, 1};
    vec nucleus7  {1, 0, 1};
    vec nucleus8  {-1, 0, 0};
    vec nucleus9  {0, -1, 0};
    vec nucleus10 {0, 0, -1};
    vec nucleus11 {-1, -1, 0};
    vec nucleus12 {0, -1, -1};
    vec nucleus13 {-1, 0, -1};
    vec nucleus14 {-1, 1, 0};
    vec nucleus15 {0, 1, -1};
    vec nucleus16 {-1, 0, 1};
    vec nucleus17 {1, -1, 0};
    vec nucleus18 {0, -1, 1};
    vec nucleus19 {1, 0, -1};


    System* system = new System(1);
    Hydrogen* Hminus = new Hydrogen("6-311++G(2d,2p)", nucleus1);
    Hminus->setNumberOfElectrons(2);
    system->addAtom(Hminus);

    /*system->addAtom(new Hydrogen("3-21++G", nucleus2));
    system->addAtom(new Hydrogen("3-21++G", nucleus2));
    system->addAtom(new Hydrogen("3-21++G", nucleus3));
    system->addAtom(new Hydrogen("3-21++G", nucleus4));
    system->addAtom(new Hydrogen("3-21++G", nucleus5));
    system->addAtom(new Hydrogen("3-21++G", nucleus6));
    system->addAtom(new Hydrogen("3-21++G", nucleus7));
    system->addAtom(new Hydrogen("3-21++G", nucleus8));
    system->addAtom(new Hydrogen("3-21++G", nucleus9));
    system->addAtom(new Hydrogen("3-21++G", nucleus10));
    system->addAtom(new Hydrogen("3-21++G", nucleus11));
    system->addAtom(new Hydrogen("3-21++G", nucleus12));
    system->addAtom(new Hydrogen("3-21++G", nucleus13));
    system->addAtom(new Hydrogen("3-21++G", nucleus14));
    system->addAtom(new Hydrogen("3-21++G", nucleus15));
    system->addAtom(new Hydrogen("3-21++G", nucleus16));
    system->addAtom(new Hydrogen("3-21++G", nucleus17));
    system->addAtom(new Hydrogen("3-21++G", nucleus18));
    system->addAtom(new Hydrogen("3-21++G", nucleus19));
    */

    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    double result = solver.solve(1e-14, 1e4);
    cout << solver.dumpBasisToFile() << endl;
    //assert(std::fabs(-1.131284349300591 - result) < 1e-15);
    return 0;
}

// -13.25093 eV  UHF H-
// -14.348   eV      H- (http://nist.gov/data/PDFfiles/jpcrd68.pdf)
// -13.6     eV      H
