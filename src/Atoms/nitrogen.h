#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Nitrogen : public Atom {
private:
    void basis_321G();
    void basis_631ppGss();

public:
    Nitrogen(std::string basisName, arma::vec position);
};

