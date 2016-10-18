#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <unordered_map>
#include <string>
#include <functional>
#include <fstream>


class IntegralTable {
protected:
    struct Key {
        int i,j,k,l;
        bool operator == ( const Key& other) const {
            return (i==other.i) && (j==other.j) && (k==other.k) && (l==other.l);
        }
    };


    struct KeyHash {
        size_t operator()(const Key& ijkl) const {
            return std::hash<int>()(ijkl.i) ^ std::hash<int>()(ijkl.j) ^
                   std::hash<int>()(ijkl.k) ^ std::hash<int>()(ijkl.l);
        }
    };

private:
    std::unordered_map<Key,double,KeyHash> m_hashMap;

public:
    IntegralTable();
    double getIntegral  (int i, int j, int k, int l);
    void   inputIntegral(int i, int j, int k, int l, double integral);
    void   printAllIntegrals();
    void   printTableToFile(std::string fileName);
    bool   readTableFromFile(std::string fileName);
    void   createTwoBodyTable(std::string   fileName,
                              int           type,
                              int           numberOfBasisFunctions,
                              int           numberOfIntegrationPoints=(int)1e7);
};
