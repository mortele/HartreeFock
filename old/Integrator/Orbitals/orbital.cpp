#include "orbital.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

Orbital::Orbital(int dimensions, int numberOfQuantumNumbers) {
    m_dimensions             = dimensions;
    m_numberOfQuantumNumbers = numberOfQuantumNumbers;
}

int* Orbital::mapToOrbitals(int p, int type) {
    int* quantumNumbers = nullptr;

    // Harmonic oscillator 2D orbitals
    if (type==0) {
        quantumNumbers = new int[2];
        switch (p) {
            case 0:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = 0;
                break;
            case 1:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = -1;
                break;
            case 2:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = 1;
                break;
            case 3:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = -2;
                break;
            case 4:
                quantumNumbers[0] = 1;
                quantumNumbers[1] = 0;
                break;
            case 5:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = 2;
                break;
            case 6:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = -3;
                break;
            case 7:
                quantumNumbers[0] = 1;
                quantumNumbers[1] = -1;
                break;
            case 8:
               quantumNumbers[0] = 1;
               quantumNumbers[1] = 1;
                break;
            case 9:
                quantumNumbers[0] = 0;
                quantumNumbers[1] = 3;
                break;
            default:
                cout << "Invalid orbital <" << p << ">." << endl;
                break;
        }

    // Hydrogen 3D orbitals
    } else if (type==1) {
        quantumNumbers = new int[3];
        switch (p) {
            case 0:
                quantumNumbers[0] = 1;
                quantumNumbers[1] = 0;
                quantumNumbers[2] = 0;
                break;
            case 1:
                quantumNumbers[0] = 2;
                quantumNumbers[1] = 0;
                quantumNumbers[2] = 0;
                break;
            case 2:
                quantumNumbers[0] = 2;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = -1;
                break;
            case 3:
                quantumNumbers[0] = 2;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = 0;
                break;
            case 4:
                quantumNumbers[0] = 2;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = 1;
                break;
            case 5:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 0;
                quantumNumbers[2] = 0;
                break;
            case 6:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = -1;
                break;
            case 7:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = 0;
                break;
            case 8:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 1;
                quantumNumbers[2] = 1;
                break;
            case 9:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 2;
                quantumNumbers[2] = -2;
                break;
            case 10:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 2;
                quantumNumbers[2] = -1;
                break;
            case 11:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 2;
                quantumNumbers[2] = 0;
                break;
            case 12:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 2;
                quantumNumbers[2] = 1;
                break;
            case 13:
                quantumNumbers[0] = 3;
                quantumNumbers[1] = 2;
                quantumNumbers[2] = 2;
                break;
            default:
                cout << "Invalid orbital <" << p << ">." << endl;
                break;
        }
    } else {
        cout << "Unknown orbital type, unable to map to spin orbitals." << endl;
    }
    return quantumNumbers;
}


int* Orbital::generateQuantumNumbers(int* indices,
                                     int  oneBodyOrTwoBody,
                                     int  type) {
    int*    allQuantumNumbers;
    int     numberOfQuantumNumbers;

    // Harmonic oscillator 2D orbitals
    if (type==0) {
        if (oneBodyOrTwoBody==1) {
            allQuantumNumbers = new int[4];
            numberOfQuantumNumbers = 2;
        } else if (oneBodyOrTwoBody==2) {
            allQuantumNumbers = new int[8];
            numberOfQuantumNumbers = 2;
        }

    // Hydrogen 3D orbitals
    } else if (type==1) {
        if (oneBodyOrTwoBody==1) {
            allQuantumNumbers = new int[6];
            numberOfQuantumNumbers = 3;
        } else if (oneBodyOrTwoBody==2) {
            allQuantumNumbers = new int[6];
            numberOfQuantumNumbers = 3;
        }
    }

    int iMax = oneBodyOrTwoBody==1 ? 2 : 4;
    for (int i=0; i<iMax; i++) {
        int* quantumNumbers = Orbital::mapToOrbitals(indices[i], type);

        for (int j=0; j<numberOfQuantumNumbers; j++) {
            allQuantumNumbers[i*numberOfQuantumNumbers+j] = quantumNumbers[j];
        }
    }
    return allQuantumNumbers;
}

int Orbital::factorial(int n) {
    int factorial = 1;
    for (int i=1; i<n+1; i++) {
        factorial *= i;
    }
    return factorial;
}

double Orbital::associatedLaguerrePolynomial(   double x,
                                                int    n,
                                                int    m) {

    if(n == 0) {
        return 1.0;
    } else if(n == 1) {
        return 1 - x + std::fabs(m);
    } else {
        return ((2*(n-1)+1 + std::fabs(m) - x)*associatedLaguerrePolynomial(x,n-1,m) - ((n-1)+std::fabs(m))*associatedLaguerrePolynomial(x,n-2,m))/((double)n);
    }
    //return (n==0) ? 1 : (1 - x + std::fabs(m));
}

double Orbital::integrandOne(double*, int*) {
    return 1.;
}

double Orbital::integrandTwo(double*, int*) {
    return 1.;
}

void Orbital::updateCoordinateScales(int* allQuantumNumbers, int numberOfQuantumNumbers) {
}

