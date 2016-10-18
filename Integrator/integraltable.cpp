#include "integraltable.h"
#include "montecarlointegrator.h"
#include "Orbitals/orbital.h"
#include "Orbitals/hydrogen3d.h"
#include "Orbitals/harmonicoscillator2d.h"

using std::cout;
using std::endl;
using std::make_pair;

IntegralTable::IntegralTable() {

}

double IntegralTable::getIntegral(int i, int j, int k, int l) {
    return m_hashMap[Key{i,j,k,l}];
}

void IntegralTable::inputIntegral(int i, int j, int k, int l, double integral) {
    m_hashMap.insert(make_pair(Key{i,j,k,l}, integral));
}

void IntegralTable::printAllIntegrals() {
    for ( auto ii = m_hashMap.begin() ; ii != m_hashMap.end() ; ii++ ) {
        cout << ii->first.i << " " << ii->first.j << " " << ii->first.k << " "
             << ii->first.l << "  "<< ii->second << endl;
    }
}

void IntegralTable::printTableToFile(std::string fileName) {
    std::ofstream outFile;
    outFile.open(fileName, std::ios::out);
    for (auto pair = m_hashMap.begin(); pair != m_hashMap.end(); pair++) {
        outFile << pair->first.i << " " << pair->first.j << " "
                << pair->first.k << " " << pair->first.l << " "
                << std::setprecision(15)
                << pair->second  << endl;
    }
    outFile.close();
}

bool IntegralTable::readTableFromFile(std::string fileName) {
    std::ifstream inFile;
    inFile.open(fileName);
    if (inFile.good() == false) {
        cout << "Unable to open table file: " << fileName << endl;
        return false;
    }
    int     i, j, k, l;
    double  integral;
    while (inFile >> i >> j >> k >> l >> integral) {
        inputIntegral(i, j, k, l, integral);
    }
    inFile.close();
    cout << "Loaded file from: " << fileName << endl;
    return true;
}

void IntegralTable::createTwoBodyTable(std::string  fileName,
                                       int          type,
                                       int          numberOfBasisFunctions,
                                       int          numberOfIntegrationPoints,
                                       MonteCarloIntegrator integrator) {

    std::string          tableFileName              = fileName;
    IntegralTable        table;

    //table.readTableFromFile("../Integrator/IntegralTables/coulomb2.dat");
    table.readTableFromFile("../IntegralTables/coulomb2.dat");
    clock_t startTime = clock();

    for (int p=0; p<numberOfBasisFunctions; p++) {
    for (int q=0; q<numberOfBasisFunctions; q++) {
    for (int r=0; r<numberOfBasisFunctions; r++) {
    for (int s=0; s<numberOfBasisFunctions; s++) {
        if (p>5 || q>5 || r>5 || s>5) {
        int quantumNumbers [] = {p,q,r,s};
        int* allQuantumNumbers = Orbital::generateQuantumNumbers(quantumNumbers, 2, type);
        int m1 = allQuantumNumbers[1];
        int m2 = allQuantumNumbers[3];
        int m3 = allQuantumNumbers[5];
        int m4 = allQuantumNumbers[7];
        if ((m1+m2)==(m3+m4)) {
            clock_t integralStart = clock();
            double I = integrator.integrateTwo(allQuantumNumbers, numberOfIntegrationPoints);
            clock_t integralFinish = clock();
            double integralTime = (integralFinish-integralStart) / ((double) CLOCKS_PER_SEC);
            double elapsedTime  = (integralFinish-startTime)     / ((double) CLOCKS_PER_SEC);
            table.inputIntegral(p,q,r,s,I);
            cout << "(p,q,r,s): " << p << ", " << q << ", " << r << ", " << s
                 << ": " << I << "  integral time: " << integralTime
                 << "  elapsed time: " << elapsedTime << endl;
            }
        }}}
        table.printTableToFile(tableFileName);
        cout << "Table dumped to file: " << tableFileName << endl;
    }}
}
