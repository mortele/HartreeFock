#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Helium : public Atom {
private:
    void basis_321G();
    void basis_6311pGss();
    void basis_6311G2df2pd();

public:
    Helium(std::string basisName, arma::vec position);
};

