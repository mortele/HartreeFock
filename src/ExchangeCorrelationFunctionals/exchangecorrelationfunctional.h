#pragma once

class ExchangeCorrelationFunctional {
protected:
    class System* m_system;

public:
    ExchangeCorrelationFunctional(class System* system);
    virtual double evaluate(double,double,double,int,int) = 0;
};

