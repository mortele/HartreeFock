#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Boron : public Atom {
private:
    void basis_321G();
    void basis_6311ppGss();

public:
    Boron(std::string basisName, arma::vec position);

};
