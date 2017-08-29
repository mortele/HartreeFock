#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>

class Hydrogen : public Atom {
private:
    void basis_321G();
    void basis_6311ppGss();

    void basis_321ppG();
    void basis_631G();
    void basis_631Gss();
    void basis_631ppGss();
    void basis_6311ppG2d2p();
    void test();

public:
    Hydrogen(std::string basisName, arma::vec position);
};

