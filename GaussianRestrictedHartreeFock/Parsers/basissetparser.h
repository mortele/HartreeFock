#pragma once
#include <string>
#include <boost/regex.hpp>
#include <armadillo>
#include "Atoms/atom.h"
#include "Parsers/filenameparser.h"

class BasisSetParser {
private:
    int             m_charge;
    int             m_numberOfElectrons;
    int             m_numberOfOrbitals;
    std::string     m_fileName;
    FileNameParser  m_fileNameParser;

    int findChargeFromAtomName(std::string atomName);
    int findNumberOfContractedOrbitals();

public:
    Atom* newAtomFromBasisSetFile(std::string atom, std::string basis, arma::vec position);
};
