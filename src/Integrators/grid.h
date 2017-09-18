#pragma once
#include <armadillo>


class Grid {
private:
    class System* m_system;
    arma::mat m_points;
    arma::vec m_weights;

public:
    Grid(class System* system);

    void createSimpleOneAtomGrid(int radialPoints, int angularPoints, double maxRadius);
    double getPoints(int i, int j) { return m_points(i,j); }
    double getWeights(int i)       { return m_weights(i); }
};

