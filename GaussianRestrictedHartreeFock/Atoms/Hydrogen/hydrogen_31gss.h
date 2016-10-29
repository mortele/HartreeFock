#pragma once
#include "Atoms/atom.h"
#include <armadillo>

class Hydrogen_31Gss : public Atom {
public:
    Hydrogen_31Gss(arma::vec position);
    std::string getInfo() { return "Hydrogen : 6-31G**"; }
};

