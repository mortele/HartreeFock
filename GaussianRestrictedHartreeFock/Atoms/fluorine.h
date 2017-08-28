#pragma once
#include "Atoms/atom.h"
#include <string>
#include <armadillo>


class Fluorine : public Atom {
    void basis_321G();
    void basis_631ppGss();

public:
    Fluorine(std::string basisName, arma::vec position);
};

#endif // FLUORINE_H
