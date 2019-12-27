#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "curve.h"
#include "distance_metrics.h"
#include "grid.h"

using namespace std;

grid::grid(double delta, int dim, int k):delta(delta), dim(dim), k(k)
{
    for(int i=0; i<k; i++)
    {
        double r = (rand()/(RAND_MAX+1.0))*dim;

        t.push_back(r);
    }

    /*printf("grid t\n");
    for(int i=0; i<k; i++)
    {
        printf("<");
        for(int j=0; j<dim; j++)
        {
            printf("%f ", t[i][j]);
        }
        printf(">\n");
    }*/
}

grid::~grid()
{

}

void grid::hash_curve(const curve& p, curve* p_g) //outputs to p_g a curve hashed through the k grids
{
    p_g->id = p.id;
    
    double idx;
    int added_elements = 0;

    for(int j=0; j<this->k; j++)
    {
        for(int i=0; i<p.length; i++)
        {
            coord c;
            c.dim = p.dim;

            for(int k=0; k<c.dim; k++)
            {
                idx = (p.data[i].c[k] - this->t[j]) / this->delta;
                double temp = (round(idx)*this->delta);
                c.c.push_back(temp);
            }

        if(added_elements++ == 0)
            p_g->add_coord(c);
        else if(p_g->data.back().cmp(c) == false)
            p_g->add_coord(c);

        }
    }
}