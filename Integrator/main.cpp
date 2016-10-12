#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include "Orbitals/orbital.h"
#include "Orbitals/harmonicoscillator2d.h"
#include "montecarlointegrator.h"
#include "ran1.h"

using std::cout;
using std::endl;

double referenceEnergy(int* p) {
    //const double reference = 1.253314137; // 1111
    //const double reference = 0.939985603; // 1212
    //const double reference = 0.3133285343; // 1221
    //const double reference = 0.744155269; // 1616
    //const double reference = 0.3035370176; // 2346
    if (p[0]==1 && p[1]==1 && p[2]==1 && p[3]==1) {
        return 1.253314137;
    } else if (p[0]==1 && p[1]==2 && p[2]==1 && p[3]==2) {
        return 0.939985603;
    } else if (p[0]==1 && p[1]==2 && p[2]==2 && p[3]==1) {
        return 0.3133285343;
    } else if (p[0]==1 && p[1]==6 && p[2]==1 && p[3]==6) {
        return 0.744155269;
    } else if (p[0]==2 && p[1]==3 && p[2]==4 && p[3]==6) {
        return 0.3035370176;
    } else {
        return NAN;
    }
}

int* mapToOrbitals(int p) {
    int* quantumNumbers = new int[2];
    switch (p) {
        case 1:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 0;
            break;
        case 2:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = -1;
            break;
        case 3:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 1;
            break;
        case 4:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = -2;
            break;
        case 5:
            quantumNumbers[0] = 1;
            quantumNumbers[1] = 0;
            break;
        case 6:
            quantumNumbers[0] = 0;
            quantumNumbers[1] = 2;
            break;
        default:
            cout << "Invalid orbital <" << p << ">." << endl;
            break;
    }
    return quantumNumbers;
}

int* generateQuantumNumbersOne(int* indices) {
    int* allQuantumNumbers = new int[4];
    for (int i=0; i<2; i++) {
        int* quantumNumbers = mapToOrbitals(indices[i]);
        allQuantumNumbers[2*i+0] = quantumNumbers[0];
        allQuantumNumbers[2*i+1] = quantumNumbers[1];
    }
    return allQuantumNumbers;
}

int* generateQuantumNumbersTwo(int* indices) {
    int* allQuantumNumbers = new int[8];
    for (int i=0; i<4; i++) {
        int* quantumNumbers = mapToOrbitals(indices[i]);
        allQuantumNumbers[2*i+0] = quantumNumbers[0];
        allQuantumNumbers[2*i+1] = quantumNumbers[1];
    }
    return allQuantumNumbers;
}

int main() {
    //int indices [] = {1,1,nan,nan};
    //int* allQuantumNumbers = generateAllQuantumNumbersOne(indices);

    int indices [] = {2,3,4,6};
    int* allQuantumNumbers = generateQuantumNumbersTwo(indices);

    MonteCarloIntegrator* MCInt = new MonteCarloIntegrator();
    MCInt->setOrbital(new HarmonicOscillator2D());
    // double I = MCInt->integrateOne(allQuantumNumbers);
    double I = MCInt->integrateTwo(allQuantumNumbers, 4e7);

    cout << "Integral:  " << I << endl;
    cout << "Reference: " << referenceEnergy(indices) << endl;
    cout << "|ref-I|:   " << std::fabs(referenceEnergy(indices)-I) << endl;
    cout << "stdDev:    " << MCInt->getStandardDeviation() << endl;

    return 0;
}
