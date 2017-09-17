#pragma once

class MonteCarloIntegrator {
private:
    long            m_seed                      = -1;
    int             m_dimensions                = 0;
    int             m_numberOfQuantumNumbers    = 0;
    bool            m_orbitalSet                = false;
    double          m_integral                  = 0;
    double          m_variance                  = 0;
    double          m_standardDeviation         = 0;
    double*         m_coordinateScales          = nullptr;
    class Orbital*  m_orbital                   = nullptr;

    void pickNewRandomCoordinates(double* singleParticleCoordinates);
    void pickAllNewRandomCoordinates(double* allCoordinates,
                                     int     numberOfParticles);

public:
    void   setOrbital       (class Orbital* orbital);
    void   setSeed          (long seed);
    double integrateOne     (int* allQuantumNumbers,
                             int  integrationPoints=(int)1e7);
    double integrateTwo     (int* allQuantumNumbers,
                             int  integrationPoints=(int)1e7);
    void   updateCoordinateScales(int* allQuantumNumbers,
                                  int  numberOfQuantumNumbers);

    double getVariance()          { return m_variance; }
    double getStandardDeviation() { return m_standardDeviation; }
};

