#pragma once
#include "Atoms/atom.h"
#include <string>

class Oxygen : public Atom {
private:
    void basis_321G();

public:
    Oxygen(std::string basisName, arma::vec position);
};

