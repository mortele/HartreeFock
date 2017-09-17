#include "numericalintegrator.h"
#include "system.h"
#include "Integrators/grid.h"

NumericalIntegrator::NumericalIntegrator(System* system) {
    m_system = system;
    m_grid   = new Grid(system);
    m_grid->createSimpleOneAtomGrid();
}
