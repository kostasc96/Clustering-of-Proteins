#include <vector>

#include "initialization.h"
#include "assignment.h"
#include "curve.h"
#include "distance_metrics.h"
#include "hash.h"

void lloyd_assignment(const vector<curve*>& dataset, const vector<curve*>& centroids, vector<int>& clusters, double& obj_function)
{
	clusters.resize(dataset.size());

	obj_function = 0.0;

	double temp_dist, best_dist;
	int best_cluster;

	for(int i=0; i<dataset.size(); i++)
	{
		best_cluster = 0;
		best_dist = calc_dist_update_cache(*dataset[i], i, *centroids[0], centroids[0]->id);

		for(int j=1; j<centroids.size(); j++)
		{
			temp_dist = calc_dist_update_cache(*dataset[i], i, *centroids[j], centroids[j]->id);

			if(temp_dist < best_dist)
			{
				best_dist = temp_dist;
				best_cluster = j;
			}
		}

		obj_function += best_dist;
		clusters[i] = best_cluster;
	}
}

void reverse_assignment(const vector<curve*>& dataset, const vector<vector_curve*>& dataset_vectors, const vector<curve*>& centroids, vector<int>& clusters, double& obj_function)
{
	vector<int> conflicts;
	vector<int> flagged(dataset.size(), 0);

	clusters.resize(dataset.size());

	for(int i=0; i<clusters.size(); i++)
		clusters[i] = -1;

	obj_function = 0.0;

	double temp_dist, temp_dist2, best_dist, R = -1.0;
	int best_cluster;
	bool res = true;

	/*min distance between centers*/
	for(int i=0; i<centroids.size()-1; i++)
	{
		best_dist = calc_dist_update_cache(*centroids[i], centroids[i]->id, *centroids[i+1], centroids[i+1]->id);

		for(int j=i+2; j<centroids.size(); j++)
		{
			temp_dist = calc_dist_update_cache(*centroids[i], centroids[i]->id, *centroids[j], centroids[j]->id);

			if(temp_dist < best_dist)
				best_dist = temp_dist;
		}

		if(R == -1.0)
			R = best_dist;
		else if(R < best_dist)
			R = best_dist;
	}

	R /= 2;
	/*min distance between centers*/

	while(res == true)
	{
		res = false;

		for(int i=0; i<centroids.size(); i++)
		{
        	for(int l=0; l<L; l++) //search for each centroid in all L tables
        	{
        		bool local_res;

        		if(centroids[i]->id != -1)
            		local_res = table[l]->find_nn_r(dataset_vectors[centroids[i]->id+(l*dataset.size())], R, i, clusters, flagged, conflicts, obj_function);
            	else
            	{
            		curve temp_g(-1, centroids[i]->dim);

            		g[l]->hash_curve(*centroids[i], &temp_g);

            		vector_curve* centroid_vc = new vector_curve;
            		centroid_vc->curve_to_vector(temp_g);
            		centroid_vc->src = centroids[i];

            		local_res = table[l]->find_nn_r(centroid_vc, R, i, clusters, flagged, conflicts, obj_function);

            		delete centroid_vc;
            	}

            	if(local_res == true)
            		res = true;
        	}
		}

		for(int i=0; i<conflicts.size(); i++)
		{
			best_cluster = 0;
			best_dist = calc_dist_update_cache(*dataset[conflicts[i]], conflicts[i], *centroids[0], centroids[0]->id);

			if(clusters[conflicts[i]] == 0)
				obj_function -= best_dist;

			for(int j=1; j<centroids.size(); j++)
			{
				temp_dist = calc_dist_update_cache(*dataset[conflicts[i]], conflicts[i], *centroids[j], centroids[j]->id);

				if(clusters[conflicts[i]] == j)
					obj_function -= temp_dist;

				if(temp_dist < best_dist)
				{
					best_dist = temp_dist;
					best_cluster = j;
				}
			}

			obj_function += best_dist;
			clusters[conflicts[i]] = best_cluster;
		}

		for(int i=0; i<dataset.size(); i++)
		{
			if(flagged[i] == 1)
				flagged[i] = 2;
		}

		conflicts.clear();

		R *= 2;
	}

	for(int i=0; i<dataset.size(); i++)
	{
		if(flagged[i] != 2)
		{
			best_cluster = 0;
			best_dist = calc_dist_update_cache(*dataset[i], i, *centroids[0], centroids[0]->id);

			for(int j=1; j<centroids.size(); j++)
			{
				temp_dist = calc_dist_update_cache(*dataset[i], i, *centroids[j], centroids[j]->id);

				if(temp_dist < best_dist)
				{
					best_dist = temp_dist;
					best_cluster = j;
				}
			}

			obj_function += best_dist;
			clusters[i] = best_cluster;
		}
	}
}