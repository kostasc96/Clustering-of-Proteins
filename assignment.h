#ifndef _ASSIGN_H
#define _ASSIGN_H

#include <vector>

#include "curve.h"
#include "grid.h"
#include "hash.h"

extern vector<grid*> g;
extern vector<LSH_table*> table;

void lloyd_assignment(const vector<curve*>&, const vector<curve*>&, vector<int>&, double&);
void reverse_assignment(const vector<curve*>&, const vector<vector_curve*>&, const vector<curve*>&, vector<int>&, double&);

#endif // _ASSIGN_H