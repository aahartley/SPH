#include "GridSearch.h"

GridSearch::GridSearch(int width, float radius)
{
    this->radius = radius;
    dimension = std::round(((width/SCALE)/radius));
    size = dimension * dimension;
    //cells = std::vector< std::vector<float> > (size);
}
void GridSearch::InitGrid(std::vector<Particle>& particles)
{
    cells = std::vector< std::vector<float> > (size);
    int count = 0;
    //std::cout << "Test\n";
    for(int j = 0; j < dimension; j++) //row
    {
        for(int i = 0; i < dimension; i++) //col
        {
            float top = (j*radius) + radius;
            float bot = (j*radius);
            float left = i * radius;
            float right = (i*radius) + radius;
            for(int k = 0; k < (int)particles.size(); k++)
            {
                if(particles[k].GetPosN().X() <= right && particles[k].GetPosN().X() >= left
                    && particles[k].GetPosN().Y() <= top && particles[k].GetPosN().Y() >= bot)
                    {
                        cells[j*dimension + i].push_back(k);
                        count++;
                        //std::cout << count << '\n';
                    }
            }
        }
    }

    neighbors = std::vector< std::vector<float> > (particles.size());

    for(int j = 0; j < dimension; j++) //rows
    {
        for(int i = 0; i < dimension; i++) //cols
        {
            for(int k = 0; k < (int)cells[j*dimension + i].size(); k++) // every particle in cell find neighbors
            {
                int index = cells[j*dimension + i][k];
                std::set<float> unique;
                //std::cout << index << '\n';
                if(j == dimension-1)
                {
                    if(i == dimension-1)
                    {
                        for(int m = 0; m < (int)cells[j*dimension + i].size(); m++)
                        {
                            //neighbors[index].push_back(cells[j*dimension + i][m]);
                            unique.insert(cells[j*dimension + i][m]);
                        }
                        for(int m = 0; m < (int)cells[(j-1)*dimension + (i)].size(); m++)
                        {
                            //neighbors[index].push_back(cells[(j-1)*dimension + i][m]);
                            unique.insert(cells[(j-1)*dimension + i][m]);

                        }
                        for(int m = 0; m < (int)cells[(j-1)*dimension + (i-1)].size(); m++)
                        {
                           //neighbors[index].push_back(cells[(j-1)*dimension + (i-1)][m]);
                            unique.insert(cells[(j-1)*dimension + (i-1)][m]);

                        }
                        for(int m = 0; m < (int)cells[(j)*dimension + (i-1)].size(); m++)
                        {
                            unique.insert(cells[j*dimension + (i-1)][m]);
                            //neighbors[index].push_back(cells[j*dimension + (i-1)][m]);
                        }
                    }
                    else
                    {
                        for(int m = 0; m < (int)cells[j*dimension + i].size(); m++)
                        {
                            //neighbors[index].push_back(cells[j*dimension + i][m]);
                            unique.insert(cells[j*dimension + i][m]);

                        }
                        for(int m = 0; m < (int)cells[(j-1)*dimension + i].size(); m++)
                        {
                            //neighbors[index].push_back(cells[(j-1)*dimension + i][m]);
                            unique.insert(cells[(j-1)*dimension + i][m]);

                        }
                        for(int m = 0; m < (int)cells[(j-1)*dimension + (i+1)].size(); m++)
                        {
                            //neighbors[index].push_back(cells[(j-1)*dimension + (i+1)][m]);
                            unique.insert(cells[(j-1)*dimension + (i+1)][m]);

                        }
                        for(int m = 0; m < (int)cells[j*dimension + (i+1)].size(); m++)
                        {
                            //neighbors[index].push_back(cells[j*dimension + (i+1)][m]);
                            unique.insert(cells[j*dimension + (i+1)][m]);

                        }
               
                    }

                }
                else if (i == dimension-1)
                {
                    for(int m = 0; m < (int)cells[j*dimension + i].size(); m++)
                    {
                        //neighbors[index].push_back(cells[j*dimension + i][m]);
                        unique.insert(cells[j*dimension + i][m]);

                    }
                    for(int m = 0; m < (int)cells[j*dimension + (i-1)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[j*dimension + (i-1)][m]);
                        unique.insert(cells[j*dimension + (i-1)][m]);

                    }
                    for(int m = 0; m < (int)cells[(j+1)*dimension + (i-1)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[(j+1)*dimension + (i-1)][m]);
                        unique.insert(cells[(j+1)*dimension + (i-1)][m]);

                    }
                    for(int m = 0; m < (int)cells[(j+1)*dimension + (i)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[(j+1)*dimension + (i)][m]);
                        unique.insert(cells[(j+1)*dimension + (i)][m]);

                    }                   
                }
                else
                {
                    for(int m = 0; m < (int)cells[j*dimension + i].size(); m++)
                    {
                        //neighbors[index].push_back(cells[j*dimension + i][m]);
                        unique.insert(cells[j*dimension + i][m]);

                    }
                    for(int m = 0; m < (int)cells[j*dimension + (i+1)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[j*dimension + (i+1)][m]);
                        unique.insert(cells[j*dimension + (i+1)][m]);

                    }
                    for(int m = 0; m < (int)cells[(j+1)*dimension + (i+1)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[(j+1)*dimension + (i+1)][m]);
                        unique.insert(cells[(j+1)*dimension + (i+1)][m]);

                    }
                    for(int m = 0; m < (int)cells[(j+1)*dimension + (i)].size(); m++)
                    {
                        //neighbors[index].push_back(cells[(j+1)*dimension + (i)][m]);
                        unique.insert(cells[(j+1)*dimension + (i)][m]);

                    }                
                }
                  
                neighbors[index].clear();         
                neighbors[index].assign(unique.begin(), unique.end());
                //std::cout << index << ": ";
                //for(int z = 0; z < neighbors[index].size(); z++)
                    //std::cout << neighbors[index][z] << ' ';
                //std::cout << '\n';
            }
        }
    }
}
