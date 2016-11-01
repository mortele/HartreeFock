#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Hydrogen : public Atom {
private:
    void basis_321G();
    void basis_321Gpp();
    void basis_631G();
    void basis_631Gss();

public:
    Hydrogen(std::string basisName, arma::vec position);
};

