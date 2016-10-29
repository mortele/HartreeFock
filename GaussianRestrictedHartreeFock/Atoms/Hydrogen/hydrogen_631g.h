#pragma once
#include "Atoms/atom.h"
#include <armadillo>
#include <vector>
#include <string>

class Hydrogen_631G : public Atom {
public:
    Hydrogen_631G(arma::vec position);
    std::string getInfo() { return "Hydrogen : 6-31G"; }
};

