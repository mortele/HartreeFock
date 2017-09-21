#pragma once

class ExchangeCorrelationFunctional {
protected:
    class System* m_system;

public:
    ExchangeCorrelationFunctional(class System* system);
    virtual double evaluateEnergy(double,double,double,int,int) = 0;
    virtual double evaluatePotential(double,double,double,int,int) = 0;
};

