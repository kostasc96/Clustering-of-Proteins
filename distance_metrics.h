#ifndef _DIST_METRICS_H
#define _DIST_METRICS_H

#include <vector>
#include "curve.h"

using namespace std;

extern vector< vector<double> > cached_distances;

double discrete_frechet_distance(const curve&, const curve&, curve&, const bool&);
double dynamic_time_warping(const curve&, const curve&);

double euclid_norm(const coord&, const coord&);

double calc_dist_update_cache(const curve&, const int&, const curve&, const int&);

double c_RMSD(const curve&, const curve&);

#endif // _DIST_METRICS_H
