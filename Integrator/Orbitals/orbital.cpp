#include "orbital.h"

Orbital::Orbital(int dimensions, int numberOfQuantumNumbers) {
    m_dimensions             = dimensions;
    m_numberOfQuantumNumbers = numberOfQuantumNumbers;
}

int* Orbital::mapQuantumNumbers(int singleQuantumNumber) {
    int* quantumNumbers = new int[1];
    quantumNumbers[0] = singleQuantumNumber;
    return quantumNumbers;
}

double Orbital::integrandOne(double* allCoordinates,
                             int*    allQuantumNumbers) {
    return 1.;
}

double Orbital::integrandTwo(double* allCoordinates,
                             int*    allQuantumNumbers) {
    return 1.;
}
