#include "tester.h"
#include <iostream>

using std::cout;
using std::endl;

Tester::Tester() {
    m_integralTester = new IntegralTester();
}

bool Tester::runAllTests() {
    cout << "┏━━━━━━━━━━━━┓" << endl;
    cout << "┃ Running  all tests ┃" << endl;
    cout << "┗━━━━━━━━━━━━┛" << endl;
    bool integralTests = runAllIntegralTests();
}

bool Tester::runAllIntegralTests() {
    cout << "Running: Integral tests" << endl;
    return m_integralTester->runAllTests();
}
