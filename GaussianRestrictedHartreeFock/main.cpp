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
#include "Tests/tester.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {
    //Examples::Hm();
    //Examples::He();
    //Examples::HeHp();
    //Examples::H2();
    Examples::H20();

    Tester tests;
    bool test = tests.runAllTests(true);
    if (!test) cout << "Some tests FAILED." << endl;
    return 0;
}
