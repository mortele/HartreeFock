#pragma once
#include "integraltester.h"

class Tester {
private:
    IntegralTester* m_integralTester;

public:
    Tester();

    bool runAllTests(bool silent=false);
    bool runAllIntegralTests(bool silent=false);
};
