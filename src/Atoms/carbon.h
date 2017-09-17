#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>


class Carbon : public Atom {
private:
    void basis_321G();
    void basis_6311ppGss();

public:
    Carbon(std::string basisName, arma::vec position);
};

