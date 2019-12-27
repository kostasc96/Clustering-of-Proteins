#include <vector>

#include "update.h"
#include "curve.h"
#include "distance_metrics.h"
#include "tree.h"

void compute_objective_function(const vector<curve*>&, const vector<int>&, const vector<curve*>&, double&);
void compute_second_cluster(const vector<curve*>&, const vector<int>&, const vector<curve*>&, vector<int>&);

bool pam_update(const vector<curve*>& dataset, const vector<int>& clusters, double& obj_function, vector<curve*>& centroids)
{
	bool ret = false;

	double dist_nc;
	double dist_sc;
	double dist_oc;

	double local_obj_function = obj_function;

	vector<int> second_cluster;
	compute_second_cluster(dataset, clusters, centroids, second_cluster);

	double final_obj = -1.0;
	int final_centroid;
	int for_cluster;

	for(int i=0; i<centroids.size(); i++)
	{
		double best_obj = -1.0;
		int best_centroid;

		vector<curve*> local_centroids = centroids;

		for(int j=0; j<clusters.size(); j++)
		{
			if(j == centroids[i]->id || clusters[j] != i)
				continue;

			vector<int> local_clusters = clusters;

			double local_obj = -1.0;

			for(int k=0; k<clusters.size(); k++)
			{
				dist_oc = calc_dist_update_cache(*dataset[k], k, *centroids[clusters[k]], centroids[clusters[k]]->id);

				if(clusters[k] == i)
				{
					dist_nc = calc_dist_update_cache(*dataset[k], k, *dataset[j], j);
					dist_sc = calc_dist_update_cache(*dataset[k], k, *centroids[second_cluster[k]], centroids[second_cluster[k]]->id);

					if(dist_nc > dist_sc)
						local_clusters[k] = second_cluster[k];
				}
				else
				{
					dist_nc = calc_dist_update_cache(*dataset[k], k, *dataset[j], j);

					if(dist_nc < dist_oc)
						local_clusters[k] = clusters[j];
				}
			}

			local_centroids[i] = dataset[j];
			if(local_centroids[i]->id != centroids[i]->id)
				compute_objective_function(dataset, local_clusters, local_centroids, local_obj);

			if(best_obj == -1.0)
			{
				best_obj = local_obj;
				best_centroid = j;
			}
			else if(local_obj < best_obj)
			{
				best_obj = local_obj;
				best_centroid = j;
			}
		}

		if(final_obj == -1.0 && best_obj != -1.0)
		{
			final_obj = best_obj;
			final_centroid = best_centroid;
			for_cluster = i;
		}
		else if(best_obj != -1.0 && best_obj < final_obj)
		{
			final_obj = best_obj;
			final_centroid = best_centroid;
			for_cluster = i;			
		}
	}

	if(final_obj != -1.0 && final_obj < local_obj_function)
	{
		centroids[for_cluster] = dataset[final_centroid];
		obj_function = final_obj;
		ret = true;
	}

	return ret;
}

void compute_objective_function(const vector<curve*>& dataset, const vector<int>& clusters, const vector<curve*>& centroids, double& obj_function)
{
	obj_function = 0.0;

	double dist;

	for(int i=0; i<clusters.size(); i++)
		obj_function += calc_dist_update_cache(*dataset[i], dataset[i]->id, *centroids[clusters[i]], centroids[clusters[i]]->id);
}

void compute_second_cluster(const vector<curve*>& dataset, const vector<int>& clusters, const vector<curve*>& centroids, vector<int>& second_cluster)
{
	second_cluster.resize(dataset.size());

	double best_dist, temp_dist;

	for(int i=0; i<clusters.size(); i++)
	{
		best_dist = -1.0;

		for(int j=0; j<centroids.size(); j++)
		{
			if(j == clusters[i])
				continue;
			
			temp_dist = calc_dist_update_cache(*dataset[i], i, *centroids[j], centroids[j]->id);

			if(best_dist == -1.0)
			{
				second_cluster[i] = j;
				best_dist = temp_dist;
			}
			else if(temp_dist < best_dist)
			{
				second_cluster[i] = j;
				best_dist = temp_dist;
			}					
		}
	}
}

bool mean_update(const vector<curve*>& dataset, const vector<int>& clusters, vector<curve*>& centroids)
{
    vector<curve*> in_cluster;

    bool ret = false;

    for(int i=0; i<centroids.size(); i++)
	{
		in_cluster.clear();

		for(int j=0; j<clusters.size(); j++)
			if(clusters[j] == i)
				in_cluster.push_back(dataset[j]);

		curve* prev_centroid = centroids[i];

		centroids[i] = calculate_mean(in_cluster);

		if(prev_centroid->cmp(*centroids[i]) == false)
			ret = true;

		if(prev_centroid->id == -1)
			delete prev_centroid;
	}

	return ret;
}