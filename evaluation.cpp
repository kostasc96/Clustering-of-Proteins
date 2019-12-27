#include "evaluation.h"

#include <vector>
#include <cstdio>

#include "curve.h"
#include "distance_metrics.h"
#include "update.h"

void silhouette(const vector<curve*>& dataset, const vector<int>& clusters, const vector<curve*>& centroids, vector<double>& s)
{
	vector<int> second_cluster(dataset.size());
	compute_second_cluster(dataset, clusters, centroids, second_cluster);

	s.resize(centroids.size(), 0.0);

	//printf("s: ");
	//for(int i=0; i<s.size(); i++)
	//	printf("%f ", s[i]);
	//printf("\n");
	//getchar();

	vector<int> n(centroids.size(), 0);

	for(int i=0; i<dataset.size(); i++)
	{
		double a = 0.0, b = 0.0;

		for(int j=0; j<clusters.size(); j++)
		{
			if(i == j)
				continue;

			if(clusters[i] == clusters[j])
				a += calc_dist_update_cache(*dataset[i], i, *dataset[j], j);
			else if(second_cluster[i] == clusters[j])
				b += calc_dist_update_cache(*dataset[i], i, *dataset[j], j);
		}

		s[clusters[i]] += (a >= b) ? (b-a)/a : (b-a)/b;
		n[clusters[i]]++;

		//printf("s(%d): %f, s[%d]: %f, n[%d]: %d second_cluster %d\n", i, (a >= b) ? (b-a)/a : (b-a)/b, clusters[i], s[clusters[i]], clusters[i], n[clusters[i]], second_cluster[i]);
	}

	for(int i=0; i<s.size(); i++)
	{
		s[i] /= n[i];
		//printf("s after div %f\n", s[i]);
	}
}