#ifndef GRIDSEARCH_H
#define GRIDSEARCH_H
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "Particle.h"
class GridSearch
{
public:
    GridSearch(){};
    GridSearch(int width, float radius);
    void InitGrid(std::vector<Particle>& particles);

    std::vector< std::vector<float> >& GetNeighbors(){ return neighbors;}

private:
    float radius;
    float diameter;
    int dimension;
    int size;
    std::vector< std::vector<float> > cells;
    std::vector< std::vector<float> > neighbors;
};


#endif