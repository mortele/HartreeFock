#include "tester.h"
#include <iostream>

using std::cout;
using std::endl;

Tester::Tester() {
    m_integralTester = new IntegralTester();
}

bool Tester::runAllTests(bool silent) {
    if (!silent) cout << "┏━━━━━━━━━━━━┓" << endl;
    if (!silent) cout << "┃ Running  all tests ┃" << endl;
    if (!silent) cout << "┗━━━━━━━━━━━━┛" << endl;
    bool integralTests = runAllIntegralTests(silent);
}

bool Tester::runAllIntegralTests(bool silent) {
    if (!silent) cout << "Running: Integral tests" << endl;
    return m_integralTester->runAllTests(silent);
}
