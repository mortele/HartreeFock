#include "montecarlointegrator.h"
#include <cmath>
#include "Orbitals/orbital.h"
#include "ran1.h"

void MonteCarloIntegrator::pickNewRandomCoordinates(double* singleParticleCoordinates) {
    for (int i = 0; i < m_dimensions; i++) {
        singleParticleCoordinates[i] = ran1(&m_seed) * m_coordinateScales[i];
    }
}

void MonteCarloIntegrator::pickAllNewRandomCoordinates(double*  allCoordinates,
                                                       int      numberOfParticles) {
    double singleParticleCoordinates[m_dimensions];
    for (int i=0; i<numberOfParticles; i++) {
        pickNewRandomCoordinates(singleParticleCoordinates);
        for (int j=0; j<m_dimensions; j++) {
            allCoordinates[i*m_dimensions + j] = singleParticleCoordinates[j];
        }
    }
}

void MonteCarloIntegrator::setOrbital(Orbital* orbital) {
    m_orbital                = orbital;
    m_dimensions             = orbital->getNumberOfDimensions();
    m_numberOfQuantumNumbers = orbital->getNumberOfQuantumNumbers();
    m_coordinateScales       = orbital->getCoordinateScales();
    m_orbitalSet             = true;
}

void MonteCarloIntegrator::setSeed(long seed) {
    m_seed = seed;
}

double MonteCarloIntegrator::integrateTwo(int* allQuantumNumbers,
                                          int  integrationPoints) {
    int     totalNumberOfCoordinates    = 2 * m_dimensions;
    double  allCoordinates[totalNumberOfCoordinates];

    m_integral = 0;
    m_variance = 0;
    for (int i=0; i<integrationPoints; i++) {
        pickAllNewRandomCoordinates(allCoordinates, 2);
        double term = m_orbital->integrandTwo(allCoordinates, allQuantumNumbers);
        m_integral += term;
        m_variance += term * term;
    }
    for (int i=0; i<2; i++) {
        for (int j=0; j<m_dimensions; j++) {
            m_integral *= m_coordinateScales[j];
            m_variance *= m_coordinateScales[j] * m_coordinateScales[j];
        }
    }
    m_integral /= ((double) integrationPoints);
    m_variance /= ((double) integrationPoints * integrationPoints);
    m_standardDeviation = std::sqrt(m_variance);
    return m_integral;
}

double MonteCarloIntegrator::integrateOne(int* allQuantumNumbers,
                                          int  integrationPoints) {
    int     totalNumberOfCoordinates    = m_dimensions;
    double  allCoordinates[totalNumberOfCoordinates];

    m_integral = 0;
    m_variance = 0;
    for (int i=0; i<integrationPoints; i++) {
        pickAllNewRandomCoordinates(allCoordinates, 1);
        double term = m_orbital->integrandOne(allCoordinates, allQuantumNumbers);
        m_integral += term;
        m_variance += term * term;
    }
    for (int j=0; j<m_dimensions; j++) {
        m_integral *= m_coordinateScales[j];
        m_variance *= m_coordinateScales[j] * m_coordinateScales[j];
    }
    m_integral /= ((double) integrationPoints);
    m_variance /= ((double) integrationPoints * integrationPoints);
    m_standardDeviation = std::sqrt(m_variance);
    return m_integral;
}
