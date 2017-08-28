#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Neon : public Atom {
private:
    void basis_321G();
    void basis_631ppGss();

public:
    Neon(std::string basisName, arma::vec position);
};

