#pragma once

class ExchangeCorrelationFunctional {
protected:
    class System* m_system;

public:
    ExchangeCorrelationFunctional(class System* system);
    virtual double evaluateEnergy(double) = 0;
    virtual double evaluatePotential(double) = 0;
};

