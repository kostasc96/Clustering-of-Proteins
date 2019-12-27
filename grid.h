#ifndef GRID_H_
#define GRID_H_

#include <vector>

#include "distance_metrics.h"
#include "curve.h"

using namespace std;

class grid
{
    public:
        grid(double, int, int);
        ~grid();

        void hash_curve(const curve&, curve*);

    private:
        double delta;
        int dim;
        int k;
        vector<double> t;
};

#endif // GRID_H_