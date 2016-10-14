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
    return true;
}
