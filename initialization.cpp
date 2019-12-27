#include <cstdio>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "curve.h"
#include "distance_metrics.h"
#include "initialization.h"

using namespace std;

int binary_search(const vector<double>&, const int&, const int&, const double& x);

void random_initialization(const vector<curve*>& dataset, const int& k, vector<curve*>& centroids)
{
	centroids.resize(k, NULL);

	int rnd;

	bool found;

	for(int i=0; i<k; i++)
	{
		do
		{
			found = false;

			rnd = (rand() / (RAND_MAX + 1.0))*((dataset.size()-1)+1);

			for(int i=0; i<centroids.size(); i++)
				if(centroids[i] != NULL && centroids[i]->id == rnd)
					found = true;

		}while(found == true);

		centroids[i] = dataset[rnd];
	}
}

void plusplus_initialization(const vector<curve*>& dataset, const int& k, vector<curve*>& centroids)
{
	centroids.resize(k, NULL);
	int t = (rand() / (RAND_MAX + 1.0))*((dataset.size()-1)+1);

	centroids[0] = dataset[t];

	vector< vector<double> > min_dists(k);

	vector<bool> if_centroid(dataset.size(), false);
	if_centroid[t] = true;

	for(int i=0; i<k; i++)
		min_dists[i].resize(dataset.size());

	for(int i=0; i<dataset.size(); i++)
	{
		if(if_centroid[i] == false)
			min_dists[0][i] = calc_dist_update_cache(*dataset[i], i, *dataset[t], t);	
	}

	vector<double> P;
	vector<int> map_P;
	P.resize(dataset.size());
	map_P.resize(dataset.size());

	for(int i=1; i<k; i++)
	{
		P[0] = 0.0;
		map_P[0] = -1;

		int j;
		for(j=0; j<dataset.size(); j++)
		{
			if(if_centroid[j] == false)
			{
				P[1] = min_dists[i-1][j]*min_dists[i-1][j];
				map_P[1] = j;
				break;
			}
		}

		int p_idx = 2;
		for(int n=j+1; n<dataset.size(); n++)
		{
			if(if_centroid[n] == false)
			{
				P[p_idx] = P[p_idx-1] + (min_dists[i-1][n]*min_dists[i-1][n]);
				map_P[p_idx] = n;
				p_idx++;
			}
		}

		double x = P[0] + (rand() / (RAND_MAX + 1.0))*(P[p_idx-1]-P[0]);

		int r = binary_search(P, 1, p_idx-1, x);
		
		centroids[i] = dataset[map_P[r]];
		if_centroid[map_P[r]] = true;

		if(i != k-1)
		{
			for(int n=0; n<dataset.size(); n++)
			{
				if(if_centroid[n] == false)
				{
					min_dists[i][n] = calc_dist_update_cache(*dataset[n], n, *dataset[map_P[r]], map_P[r]);

					if(min_dists[i-1][n] < min_dists[i][n])
						min_dists[i][n] = min_dists[i-1][n];
				}
			}
		}
	}
}

int binary_search(const vector<double>& P, const int& start, const int& end, const double& x)
{
    if(start > end)
        return -1;

    const int middle = start + ((end - start) / 2);

    if(P[middle-1] < x && x <= P[middle])
        return middle;
    else if(P[middle-1] >= x)
        return binary_search(P, start, middle-1, x);

    return binary_search(P, middle+1, end, x);
}