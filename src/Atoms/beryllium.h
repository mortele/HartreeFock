#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Beryllium : public Atom {
private:
    void basis_321G();
    void basis_631ppGss();
    void basis_STO6G();

public:
    Beryllium(std::string basisName, arma::vec position);
};

