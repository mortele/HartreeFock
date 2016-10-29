#pragma once
#include "Atoms/atom.h"
#include <armadillo>
#include <vector>
#include <string>

class Hydrogen_321G : public Atom {
public:
    Hydrogen_321G(arma::vec position);
    std::string getInfo() { return "Hydrogen : 3-21G"; }
};
