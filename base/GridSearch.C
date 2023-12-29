#include "GridSearch.h"

GridSearch::GridSearch(int width, float radius)
{
    this->radius = radius;
    this->diameter = radius*2;
    dimension = std::round(((width/SCALE)/(diameter)));
    size = dimension * dimension;
}
void GridSearch::InitGrid(std::vector<Particle>& particles)
{
    cells = std::vector< std::vector<float> > (size);
    //int count = 0;
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < dimension; j++) //row
    {
        for(int i = 0; i < dimension; i++) //col
        {
            float top = (j*diameter) + diameter;
            float bot = (j*diameter);
            float left = i * diameter;
            float right = (i*diameter) + diameter;
            for(int k = 0; k < (int)particles.size(); k++)
            {
                if(particles[k].GetPosN().X() <= right && particles[k].GetPosN().X() >= left
                    && particles[k].GetPosN().Y() <= top && particles[k].GetPosN().Y() >= bot)
                    {
                        #pragma omp critical
                        cells[j*dimension + i].push_back(k);
                        //count++;
                    }
            }
        }
    }
    //std::cout << count << '\n';


    neighbors = std::vector< std::vector<float> > (particles.size());
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < dimension; j++) //rows
    {
        for(int i = 0; i < dimension; i++) //cols
        {
            float bot = (j*diameter);
            float left = i * diameter;
            for(int k = 0; k < (int)cells[j*dimension + i].size(); k++) // every particle in cell find neighbors
            {
                int index = cells[j*dimension + i][k]; //particle index
                std::vector<int> indices(4);
                std::set<float> unique;
                float distanceY = std::abs(particles[index].GetPosN().Y() - bot);
                float distanceX = std::abs(particles[index].GetPosN().X() - left);
                int maxi, mini, maxy, miny;
                int index1, index2, index3, index4;
                if(distanceX >  radius) // go right
                {
                    maxi = i + 1;
                    mini = i;
                }
                else if (distanceX < radius) // go left
                {
                    maxi = i;
                    mini = i-1;
                }
                else //distance == radius
                {
                    maxi = i;
                    mini = i;
                }
                if(distanceY > radius) // go up
                {   
                    maxy = j + 1;
                    miny = j;
                }
                else if(distanceY < radius) // go down
                {
                    maxy = j;
                    miny = j - 1;
                }
                else // distance == radius
                {
                    maxy = j;
                    miny = j;
                }
                if(miny > dimension-1) miny = dimension-1;
                if(miny < 0) miny = 0;
                if(maxy > dimension-1) maxy = dimension-1;
                if(maxy < 0) maxy = 0;
                if(maxi > dimension-1) maxi = dimension-1;
                if(maxi < 0) maxi = 0;
                if(mini > dimension-1) mini = dimension-1;
                if(mini < 0) mini = 0;
                index1 = miny * dimension + mini;
                index2 = miny * dimension + maxi;
                index3 = maxy * dimension + mini;
                index4 = maxy * dimension + maxi;
            
                indices[0] = index1; indices[1] = index2; indices[2] = index3; indices[3] = index4;
                for(int i = 0; i < (int)indices.size(); i++)
                {
                    for(int j = 0; j < (int)cells[indices[i]].size(); j++)
                    {
                        if(cells[indices[i]][j] != index)
                            unique.insert(cells[indices[i]][j]);
                    }
                }
                #pragma omp critical
                neighbors[index].assign(unique.begin(), unique.end());
            }
        }
    }
}
