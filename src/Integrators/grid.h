#pragma once
#include <armadillo>


class Grid {
private:
    class System* m_system;
    arma::mat m_points;
    arma::vec m_weights;

public:
    Grid(class System* system);

    void createSimpleOneAtomGrid(int radialPoints=200, int angularPoints=200, double maxRadius=-1);
    void beckeGrid(int minimumAngularPoints, int maximumAngularPoints);

    arma::mat getPoints()  { return m_points; }
    arma::vec getWeights() { return m_weights; }
};

