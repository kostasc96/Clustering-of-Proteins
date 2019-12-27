#include <vector>
#include <cmath>

#include "tree.h"
#include "curve.h"
#include "distance_metrics.h"

using namespace std;

static int pos;

void build_tree(vector<bst>& tree, const int& index, int& last_index, const int& local_height, const int& height, const vector<curve*>& cluster)
{
	if(local_height < height)
	{
		tree[index].lindex = last_index+1;
		tree[index].rindex = last_index+2;
		tree[index].data = NULL;
		tree[tree[index].lindex].data = NULL;
		tree[tree[index].rindex].data = NULL;
		last_index += 2;

		build_tree(tree, tree[index].lindex, last_index, local_height+1, height, cluster);

		if(pos < cluster.size())
		{
			build_tree(tree, tree[index].rindex, last_index, local_height+1, height, cluster);

			tree[index].data = new curve(-1, cluster[0]->dim);
			discrete_frechet_distance(*tree[tree[index].lindex].data, *tree[tree[index].rindex].data, *tree[index].data, true);

			if(tree[tree[index].lindex].data->id == -1)
				delete tree[tree[index].lindex].data;
			if(tree[tree[index].rindex].data->id == -1)
				delete tree[tree[index].rindex].data;
		}
		else
			tree[index].data = tree[tree[index].lindex].data;
	}
	else if(local_height == height)
	{
		tree[index].data = cluster[pos++];
		tree[index].lindex = -1;
		tree[index].rindex = -1;
	}
}

curve* calculate_mean(const vector<curve*>& cluster)
{
	pos = 0;

	int h = ceil(log2(cluster.size()));
	int bst_size = pow(2, h+1) - 1;

	vector<bst> tree(bst_size);

	int last_index = 0;

	build_tree(tree, 0, last_index, 0, h, cluster);

	return tree[0].data;
}