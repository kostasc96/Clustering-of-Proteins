#ifndef _UPDATE_H
#define _UPDATE_H

#include <vector>

#include "curve.h"

using namespace std;

bool pam_update(const vector<curve*>&, const vector<int>&, double&, vector<curve*>&);
bool mean_update(const vector<curve*>&, const vector<int>&, vector<curve*>&);

void compute_second_cluster(const vector<curve*>&, const vector<int>&, const vector<curve*>&, vector<int>&);

#endif //_UPDATE_H