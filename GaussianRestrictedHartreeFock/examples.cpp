#include "examples.h"
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


using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;
using std::cos;
using std::sin;



void Examples::Hm() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};

    System* system = new System(1);
    Hydrogen* Hm = new Hydrogen("6-311++G(2d,2p)", nucleus1);
    Hm->setNumberOfElectrons(2);
    system->addAtom(Hm);

    //UnrestrictedHartreeFock solver(system);
    RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}

void Examples::He() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};

    System* system = new System(1);
    system->addAtom(new Helium("6-311+G**", nucleus1));

    //UnrestrictedHartreeFock solver(system);
    RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}

void Examples::HeHp() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};
    vec nucleus2  {0, 0, 1.4632};

    System* system = new System(2);
    system->addAtom(new Helium("6-311+G**", nucleus1));
    Hydrogen* Hm = new Hydrogen("6-311++G(2d,2p)", nucleus2);
    Hm->setNumberOfElectrons(0);
    system->addAtom(Hm);

    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}


void Examples::H2() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};
    vec nucleus2  {0, 0, 1.4};

    System* system = new System(2);
    system->addAtom(new Hydrogen("6-311++G**", nucleus1));
    system->addAtom(new Hydrogen("6-311++G**", nucleus2));
    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}


void Examples::H20() {
    boost::timer t;
    double r        = 1.809;
    double theta    = 1.824;

    vec nucleus1  {0,            0,             0};
    vec nucleus2  {r,            0,             0};
    vec nucleus3  {r*cos(theta), r*sin(theta),  0};

    System* system = new System(3);
    system->addAtom(new Oxygen  ("6-311++G**", nucleus1));
    system->addAtom(new Hydrogen("6-311++G**", nucleus2));
    system->addAtom(new Hydrogen("6-311++G**", nucleus3));
    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    solver.solve(1e-7, 1);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;

    /*for (ContractedGaussian* contracted: system->getBasis()) {
        cout << *contracted << endl;
    }*/
}

// -13.25093 eV  UHF H-
// -14.348   eV      H- (http://nist.gov/data/PDFfiles/jpcrd68.pdf)
// -13.6     eV      H






