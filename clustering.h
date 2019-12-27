#ifndef _CLUSTER_H
#define _CLUSTER_H

#include <vector>

#include "curve.h"

using namespace std;

extern FILE* output;
extern char* metric;

double I1A1U1_clustering(const vector<curve*>&, const int&);
double I1A1U2_clustering(const vector<curve*>&, vector<int>&, const int&);
double I1A2U1_clustering(const vector<curve*>&, const vector<vector_curve*>& , const int&);
double I1A2U2_clustering(const vector<curve*>&, const vector<vector_curve*>& , const int&);
double I2A1U1_clustering(const vector<curve*>&, const int&);
double I2A1U2_clustering(const vector<curve*>&, const int&);
double I2A2U1_clustering(const vector<curve*>&, const vector<vector_curve*>& , const int&);
double I2A2U2_clustering(const vector<curve*>&, const vector<vector_curve*>& , const int&);

void cluster_with_eval(double cluster_f(const vector<curve*>&, vector<int>&, const int&), const vector<curve*>&, int, bool);

void print_output_rmsd(const vector<int>&, int, double, double, bool);

#endif //_CLUSTER_H