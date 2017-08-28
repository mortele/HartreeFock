#pragma once
#include "atom.h"

class Lithium : public Atom {
private:
    void basis_321G();
    void basis_6311ppGss();
    
public:
    Lithium(std::string basisName, arma::vec position);
};
