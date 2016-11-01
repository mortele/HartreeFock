#pragma once
#include <iostream>
#include <string>

class FileNameParser {
private:
    std::string m_path = "../EMSL_Basis_Set_Exchange_Local/db/";

public:
    std::string findBasisFile(std::string atom, std::string basis);
};
