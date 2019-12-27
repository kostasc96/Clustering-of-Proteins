#include <vector>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "clustering.h"
#include "initialization.h"
#include "distance_metrics.h"
#include "curve.h"
#include "assignment.h"
#include "grid.h"
#include "hash.h"
#include "update.h"
#include "tree.h"
#include "evaluation.h"

#define MAX_TIMES 15

using namespace std;

void print_output(const char*, const vector<int>&, const vector<curve*>&, const clock_t&, const vector<double>&);

double I1A1U1_clustering(const vector<curve*>& dataset, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I1A1U1 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

	plusplus_initialization(dataset, k, centroids);

	for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    lloyd_assignment(dataset, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    int times = 0;
    bool change;
	do
	{
    	change = mean_update(dataset, clusters, centroids);

    	lloyd_assignment(dataset, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", obj_function);

        times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I1A1U1";
    print_output(algo, clusters, centroids, t, s);*/

    for(int i=0; i<centroids.size(); i++)
        delete centroids[i];

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I1A1U2_clustering(const vector<curve*>& dataset, vector<int>& clusters, const int& k)
{
	clock_t t = clock();
	//printf("Clustering with I1A1U2 algorithm\n");

	vector<curve*> centroids;
    double obj_function;
    vector<double> s;

    plusplus_initialization(dataset, k, centroids);

    /*for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);*/

    lloyd_assignment(dataset, centroids, clusters, obj_function);

    //printf("obj_function: %f\n", obj_function);

    /*for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }*/

    double pam_obj = obj_function;

    int times = 0;
    bool change;
	do
	{
    	change = pam_update(dataset, clusters, pam_obj, centroids);

    	/*for(int i=0; i<centroids.size(); i++)
        	printf("updated centroid: %d\n", centroids[i]->id);*/

    	lloyd_assignment(dataset, centroids, clusters, obj_function);

        //printf("updated obj_function: %f\n", pam_obj);

        times++;

   	} while(change == true && times <= MAX_TIMES);

    /*for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }*/

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I1A1U2";
    print_output(algo, clusters, centroids, t, s);*/

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I1A2U1_clustering(const vector<curve*>& dataset, const vector<vector_curve*>& dataset_vectors, const int& k)
{
	clock_t t = clock();	
	printf("Clustering with I1A2U1 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

	plusplus_initialization(dataset, k, centroids);

	for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    int times = 0;
    bool change;
	do
	{
    	change = mean_update(dataset, clusters, centroids);

    	reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", obj_function);

        times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I1A2U1";
    print_output(algo, clusters, centroids, t, s);*/

    for(int i=0; i<centroids.size(); i++)
        delete centroids[i];

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I1A2U2_clustering(const vector<curve*>& dataset, const vector<vector_curve*>& dataset_vectors, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I1A2U2 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

    plusplus_initialization(dataset, k, centroids);

    for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    double pam_obj = obj_function;

    int times = 0;
    bool change;
	do
	{
    	change = pam_update(dataset, clusters, pam_obj, centroids);

    	for(int i=0; i<centroids.size(); i++)
        	printf("updated centroid: %d\n", centroids[i]->id);

    	reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", pam_obj);

    	times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I1A2U2";
    print_output(algo, clusters, centroids, t, s);*/

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I2A1U1_clustering(const vector<curve*>& dataset, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I2A1U1 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

	random_initialization(dataset, k, centroids);

	for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    lloyd_assignment(dataset, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    int times = 0;
    bool change;
	do
	{
    	change = mean_update(dataset, clusters, centroids);

    	lloyd_assignment(dataset, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", obj_function);

        times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I2A1U1";
    print_output(algo, clusters, centroids, t, s);*/

    for(int i=0; i<centroids.size(); i++)
        delete centroids[i];

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I2A1U2_clustering(const vector<curve*>& dataset, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I2A1U2 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

    random_initialization(dataset, k, centroids);

    for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    lloyd_assignment(dataset, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    double pam_obj = obj_function;

    int times = 0;
    bool change;
	do
	{
    	change = pam_update(dataset, clusters, pam_obj, centroids);

    	for(int i=0; i<centroids.size(); i++)
        	printf("updated centroid: %d\n", centroids[i]->id);

    	lloyd_assignment(dataset, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", pam_obj);

        times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I2A1U2";
    print_output(algo, clusters, centroids, t, s);*/

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I2A2U1_clustering(const vector<curve*>& dataset, const vector<vector_curve*>& dataset_vectors, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I2A2U1 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

	random_initialization(dataset, k, centroids);

	for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    int times = 0;
    bool change;
	do
	{
    	change = mean_update(dataset, clusters, centroids);

    	reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", obj_function);

        times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I2A2U1";
    print_output(algo, clusters, centroids, t, s);*/

    for(int i=0; i<centroids.size(); i++)
        delete centroids[i];

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

double I2A2U2_clustering(const vector<curve*>& dataset, const vector<vector_curve*>& dataset_vectors, const int& k)
{
	clock_t t = clock();
	printf("Clustering with I2A2U2 algorithm\n");

	vector<curve*> centroids;
	vector<int> clusters;
    double obj_function;
    vector<double> s;

    random_initialization(dataset, k, centroids);

    for(int i=0; i<centroids.size(); i++)
        printf("centroid: %d\n", centroids[i]->id);

    reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    printf("obj_function: %f\n", obj_function);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    double pam_obj = obj_function;

    int times = 0;
    bool change;
	do
	{
    	change = pam_update(dataset, clusters, pam_obj, centroids);

    	for(int i=0; i<centroids.size(); i++)
        	printf("updated centroid: %d\n", centroids[i]->id);

    	reverse_assignment(dataset, dataset_vectors, centroids, clusters, obj_function);

    	printf("updated obj_function: %f\n", pam_obj);

    	times++;

	} while(change == true && times <= MAX_TIMES);

    for(int i=0; i<centroids.size(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                printf("%d ", j);
        printf("\n");
    }

    silhouette(dataset, clusters, centroids, s);

    t = clock() - t;

    /*char *algo = "I2A2U2";
    print_output(algo, clusters, centroids, t, s);*/

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
        silhouette_coefficient += s[i];

    return silhouette_coefficient /= s.size();
}

void print_output(const char *algo, const vector<int>& clusters, const vector<curve*>& centroids, const clock_t& t, const vector<double>& s)
{
	fprintf(output, "Algorithm: %s\n", algo);
	fprintf(output, "Metric: %s\n", metric);

	vector<int> c_size(centroids.size(), 0);

	for(int i=0; i<centroids.size(); i++)
	{
		for(int j=0; j<clusters.size(); j++)
			if(clusters[j] == i)
				c_size[i]++;

		if(centroids[i]->id == -1)
		{
			fprintf(output, "CLUSTER-%d {size: %d, centroid: \n", i, c_size[i]);
			for(int j=0; j<centroids[i]->length; j++)
			{
				fprintf(output, "point %d: ", j);
				centroids[i]->data[j].print_out();
				fprintf(output, "\n");
			}
			fprintf(output, "}\n");
		}
		else
			fprintf(output, "CLUSTER-%d {size: %d, centroid %d}\n", i, c_size[i], centroids[i]->id);
	}

	fprintf(output, "clustering_time: %f\n", ((double)t)/CLOCKS_PER_SEC);

	fprintf(output, "Silhouette: [");

	double silhouette_coefficient = 0.0;
	for(int i=0; i<s.size(); i++)
	{
		silhouette_coefficient += s[i];
		fprintf(output, "%f,", s[i]);
	}

	fprintf(output, "%f]\n", silhouette_coefficient /= s.size());
}

void print_output_rmsd(const vector<int>& clusters, int k, double s, double time, bool flag)
{
    if (flag == true) {
        fprintf(output, "%d\n%f\n%f\n", k, s, time);
    } else {
        fprintf(output, "%d\n%f\n", k, s);
    }

    for(int i=0; i<k; i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(clusters[j] == i)
                fprintf(output, "%d\t", j);

        fprintf(output, "\n");
    }
}

void cluster_with_eval(double cluster_f(const vector<curve*>&, vector<int>&, const int&), const vector<curve*>& dataset, int k, bool timed)
{
    clock_t t = clock();

    vector<int> clusters;
    vector<int> prev_clusters;

    double prev_s = 0, s = 0;

    while ((s-prev_s >= (prev_s*15/100) || fabs(prev_s-s) <= (prev_s*10/100)) || (prev_s < s)) {
        prev_clusters = clusters;
        prev_s = s;

        s = cluster_f(dataset, clusters, k++);
    }

    t = clock() - t;

    double time = (((double)t)/CLOCKS_PER_SEC)*1000;

    print_output_rmsd(prev_clusters, k-2, prev_s, time, timed);
}