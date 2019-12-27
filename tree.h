#ifndef _TREE_H
#define _TREE_H

#include <vector>

#include "curve.h"

using namespace std;

struct bst
{
	curve* data;
	int lindex;
	int rindex;
};

void build_tree(vector<bst>&, const int&, int&, const int&, const int&, const vector<curve*>&);
curve* calculate_mean(const vector<curve*>&);

#endif //_TREE_H