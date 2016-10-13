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
#include "integraltable.h"

using std::cout;
using std::endl;


int* mapToOrbitals(int p) {
    int* quantumNumbers = new int[2];
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

    IntegralTable table;
    table.inputIntegral(0,0,0,0,  1.253314137);
    table.inputIntegral(0,1,0,1,  0.939985603);
    table.inputIntegral(0,1,0,1,  0.3133285343);
    table.inputIntegral(0,5,0,5,  0.744155269);
    table.inputIntegral(1,2,3,5,  0.3035370176);
    //table.printAllIntegrals();

    //int q [] = {0,0,nan,nan};
    //int* allQuantumNumbers = generateAllQuantumNumbersOne(q);

    int q [] = {1,2,3,5};
    int* allQuantumNumbers = generateQuantumNumbersTwo(q);

    MonteCarloIntegrator* MCInt = new MonteCarloIntegrator();
    MCInt->setOrbital(new HarmonicOscillator2D());
    // double I = MCInt->integrateOne(allQuantumNumbers);
    double I    = MCInt->integrateTwo(allQuantumNumbers, 1e7);
    double ref  = table.getIntegral(q[0], q[1], q[2], q[3]);

    cout << "Integral:  " << I << endl;
    cout << "Reference: " << ((ref==0) ? NAN : ref) << endl;
    cout << "|ref-I|:   " << ((ref==0) ? NAN : std::fabs(ref-I)) << endl;
    cout << "stdDev:    " << MCInt->getStandardDeviation() << endl;

    return 0;
}

