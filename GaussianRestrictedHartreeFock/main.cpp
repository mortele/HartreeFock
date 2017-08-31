//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include "examples.h"
#include "Tests/tester.h"

using std::cout;
using std::endl;

int main(int, char**) {
    //Examples::Hm();
    //Examples::He();
    //Examples::HeHp();
    Examples::H2();
    //Examples::H20();
    //Examples::ValidationTable();
    //Examples::ValidationTableDissociation();
    //Examples::ValidationH2plus();

    Tester tests;
    bool test = tests.runAllTests(true);
    if (!test) cout << "Some tests FAILED." << endl;
    return 0;
}

