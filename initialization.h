#ifndef _INIT_H
#define _INIT_H

#include <vector>

#include "curve.h"

using namespace std;

void plusplus_initialization(const vector<curve*>&, const int&, vector<curve*>&);
void random_initialization(const vector<curve*>&, const int&, vector<curve*>&);

#endif // _INIT_H