#pragma once

class NumericalIntegrator {
private:
    class System*   m_system;
    class Grid*     m_grid;

public:
    NumericalIntegrator(class System* system);
    double testIntegral();
};

