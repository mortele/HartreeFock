#include "integraltable.h"

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
    outFile.open("../IntegralTables/" + fileName, std::ios::out);
    for (auto pair = m_hashMap.begin(); pair != m_hashMap.end(); pair++) {
        outFile << pair->first.i << " " << pair->first.j << " "
                << pair->first.k << " " << pair->first.l << " "
                << std::setprecision(15)
                << pair->second  << endl;
    }
    outFile.close();
}

void IntegralTable::readTableFromFile(std::string fileName) {
    std::ifstream inFile;
    inFile.open("../IntegralTables/" + fileName);

    int     i, j, k, l;
    double  integral;
    while (inFile >> i >> j >> k >> l >> integral) {
        inputIntegral(i, j, k, l, integral);
    }
    inFile.close();
}
