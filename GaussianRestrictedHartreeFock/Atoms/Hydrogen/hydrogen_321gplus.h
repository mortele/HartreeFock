#pragma once
#include "Atoms/atom.h"
#include <armadillo>
#include <vector>
#include <string>

class Hydrogen_321Gplus : public Atom {
public:
    Hydrogen_321Gplus(arma::vec position);
    std::string getInfo() { return "Hydrogen : 3-21++G"; }
};
