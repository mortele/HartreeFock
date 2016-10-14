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
    int                  numberOfBasisFunctions     = 6;
    int                  numberOfIntegrationPoints  = (int) 5e7;
    std::string          tableFileName              = "HO_2d_6.dat";
    IntegralTable        table;
    MonteCarloIntegrator MCInt;
    MCInt.setOrbital(new HarmonicOscillator2D());

    for (int p=0; p<numberOfBasisFunctions; p++) {
        for (int q=0; q<numberOfBasisFunctions; q++) {
            for (int r=0; r<numberOfBasisFunctions; r++) {
                for (int s=0; s<numberOfBasisFunctions; s++) {
                    int quantumNumbers [] = {p,q,r,s};
                    int* allQuantumNumbers = generateQuantumNumbersTwo(quantumNumbers);
                    double I = MCInt.integrateTwo(allQuantumNumbers, numberOfIntegrationPoints);
                    table.inputIntegral(p,q,r,s,I);
                    cout << "(p,q,r,s): " << p << ", " << q << ", " << r << ", " << s
                         << ": " << I << endl;
                }
                table.printTableToFile(tableFileName);
                cout << "Table dumped to file: " << tableFileName << endl;
            }
        }
    }

    return 0;
}
