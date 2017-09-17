#pragma once
#include "Atoms/atom.h"
#include <string>

class Oxygen : public Atom {
private:
    void basis_321G();
    void basis_631pGss();
    void basis_6311ppGss();

public:
    Oxygen(std::string basisName, arma::vec position);
};

