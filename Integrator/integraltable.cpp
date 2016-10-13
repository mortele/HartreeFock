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
    for ( auto ii = m_hashMap.begin() ; ii != m_hashMap.end() ; ii++ )
        cout << ii->first.i
             << " "
             << ii->first.j
             << " "
             << ii->first.k
             << " "
             << ii->first.l
             << "  "
             << ii->second
             << endl;
}
