#pragma once

class Orbital {
public:
    Orbital(int dimensions, int quantumNumbers);
    virtual double computeWavefunction(double* coordinates,
                                       int*    numberOfQuantumNumbers) = 0;
    virtual double* getCoordinateScales() = 0;
    int getNumberOfDimensions()     { return m_dimensions; }
    int getNumberOfQuantumNumbers() { return m_numberOfQuantumNumbers; }
    virtual double integrandOne(double* allCoordinates,
                                int*    allQuantumNumbers);
    virtual double integrandTwo(double* allCoordinates,
                                int*    allQuantumNumbers);

protected:
    int     m_dimensions              = 0;
    int     m_numberOfQuantumNumbers  = 0;

    virtual int*   mapQuantumNumbers(int singleQuantumNumber);
};
