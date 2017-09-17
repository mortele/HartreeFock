#include "basissetparser.h"
#include <cmath>
#include <fstream>

using std::string;
using std::cout;
using std::endl;
using std::getline;
using arma::vec;
using std::ifstream;


Atom* BasisSetParser::newAtomFromBasisSetFile(string atomName,
                                              string basisName,
                                              vec    position) {

    m_fileName = m_fileNameParser.findBasisFile(atomName, basisName);
    m_numberOfOrbitals = findNumberOfContractedOrbitals();
    m_charge = findChargeFromAtomName(atomName);
    m_numberOfElectrons = m_charge;
    Atom* atom = new Atom(position, 0, m_numberOfElectrons, m_charge);
    atom->setInfo(atomName + " : " + basisName);
    return atom;
}

int BasisSetParser::findNumberOfContractedOrbitals() {
    ifstream basisFile(m_fileName);
    if (! basisFile.good()) {
        cout << "Unable to open basis file " << m_fileName << endl;
    }
    for (string line; getline(basisFile, line); ) {
        cout << line << endl;
    }
}

int BasisSetParser::findChargeFromAtomName(std::string atomName) {
    if (atomName =="H") {
        return 1;
    } else if (atomName =="He") {
        return 2;
    } else if (atomName =="Li") {
        return 3;
    } else if (atomName =="Be") {
        return 4;
    } else if (atomName =="B") {
        return 5;
    } else if (atomName =="C") {
        return 6;
    } else if (atomName =="N") {
        return 7;
    } else if (atomName =="O") {
        return 8;
    } else if (atomName =="F") {
        return 9;
    } else if (atomName =="Ne") {
        return 10;
    } else if (atomName =="Na") {
        return 11;
    } else if (atomName =="Mg") {
        return 12;
    } else if (atomName =="Al") {
        return 13;
    } else if (atomName =="Si") {
        return 14;
    } else if (atomName =="P") {
        return 15;
    } else if (atomName =="S") {
        return 16;
    } else if (atomName =="Cl") {
        return 17;
    } else if (atomName =="Ar") {
        return 18;
    } else {
        return NAN;
    }
}







