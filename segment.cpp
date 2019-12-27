#include <cstdio>
#include <cstring>
#include <ctime>
#include <deque>
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <map>

#include "curve.h"
#include "clustering.h"
#include "hash.h"
#include "grid.h"
#include "update.h"
#include "evaluation.h"

#define CURVATURE_CEIL 3
#define LENGTH_CEIL 50

#define DELTA 1
#define MAX_R 100

using namespace std;

FILE *output;
FILE *files;

char *metric;

int L;
int K;
int k;

vector< vector<double> > cached_distances;
vector<int> r;

vector<curve*> dataset;
vector<vector_curve*> dataset_vectors;

int dataset_size;

vector<grid*> g;
vector<LSH_table*> table;

clock_t t;

double curvature(const deque<double>& line_set) {
	double a = sqrt((line_set[0]-line_set[2])*(line_set[0]-line_set[2])+(line_set[1]-line_set[3])*(line_set[1]-line_set[3]));
	double b = sqrt((line_set[0]-line_set[4])*(line_set[0]-line_set[4])+(line_set[1]-line_set[5])*(line_set[1]-line_set[5]));
	double c = sqrt((line_set[2]-line_set[4])*(line_set[2]-line_set[4])+(line_set[3]-line_set[5])*(line_set[3]-line_set[5]));

	double area = sqrt((a+b+c)*(b+c-a)*(a-b+c)*(a+b-c));

	return a*b*c/area;
}

map<pair<double,double>,int> nodes_count;

int find_intersections(const char *in) {
	FILE* input;
	if((input = fopen(in, "r")) == 0)
		return 1;

	int way_id, segment_id = 0;
	char type[32];

	double lat, lon;

	while (fscanf(input, "%d,%[^,],", &way_id, type) > 0) {
		char temp = ',';

		while (temp != '\n') {
			if(fscanf(input, "%lf,%lf%c" , &lat, &lon, &temp) == EOF)
				break;

			if (!nodes_count.count(make_pair(lat,lon))) {
				nodes_count[make_pair(lat,lon)] = 1;
			} else {
				nodes_count[make_pair(lat,lon)]++;
			}
		}
	}

	fclose(input);

	return 0;		
}

int parse(const char *in, const char*out) {
	FILE* input;
	if((input = fopen(in, "r")) == 0)
		return 1;

	FILE* output;
	if((output = fopen(out, "w")) == 0)
		return 2;

	int way_id, segment_id = 0;
	char type[32];

	double lat, lon, prev_curv = -1.0, curv = -1.0;

	deque<double> line_set, segment;

	while (fscanf(input, "%d,%[^,],", &way_id, type) > 0) {
		char temp = ',';

		line_set.clear();
		segment.clear();

		while (temp != '\n') {
			if(fscanf(input, "%lf,%lf%c" , &lat, &lon, &temp) == EOF)
				break;

			segment.push_back(lat);
			segment.push_back(lon);

			line_set.push_back(lat);
			line_set.push_back(lon);

			if (line_set.size() == 6) {
				prev_curv = curv;
				curv = curvature(line_set);

				if (((prev_curv != -1.0 && ((curv < prev_curv) ? curv : prev_curv) >= CURVATURE_CEIL) && 
					(nodes_count[make_pair(line_set[0],line_set[1])] == 1 && nodes_count[make_pair(line_set[2],line_set[3])] == 1)) || (segment.size()-4 >= LENGTH_CEIL)) {
					fprintf(output, "%d,%d,%d", segment_id++, way_id, (int)(segment.size()-4)/2);

					int count = 0;
					for (int i=0; i<segment.size()-4; i++) {
						fprintf(output, ",%.7lf", segment[i]);
						count++;
					}
					fprintf(output, "\n");

					for (int i=0; i<count; i++) {
						segment.pop_front();
					}
					line_set.pop_front();
					line_set.pop_front();
				} else {
					line_set.pop_front();
					line_set.pop_front();
				}
			}
		}

		if (segment.size() > 0) {
			fprintf(output, "%d,%d,%d", segment_id++, way_id, (int)segment.size()/2);

			for (int i=0; i<segment.size(); i++) {
				fprintf(output, ",%.7lf", segment[i]);
			}
			fprintf(output, "\n");
		}
	}

	fclose(input);
	fclose(output);

	return 0;
}

int parse_segments(const char* in, int& max_nodes) {
	FILE* input;
	if ((input = fopen(in, "r")) == 0) {
		return 1;
	}

	int dim = 2, segment_id, way_id, num_of_nodes;

	max_nodes = -1;

    coord c;
    c.dim = dim;
    c.c.resize(dim);

	while (fscanf(input, "%d,%d,%d", &segment_id, &way_id, &num_of_nodes) > 0) {
		curve* temp_c = new curve(dataset.size(), dim);
		temp_c->data.reserve(num_of_nodes);

		if (num_of_nodes > max_nodes) {
			max_nodes = num_of_nodes;
		}

		for (int i=0; i<num_of_nodes; i++) {
			fscanf(input, ",%lf,%lf", &c.c[0], &c.c[1]);

			temp_c->add_coord(c);
		}
		fscanf(input, "\n");

		dataset.push_back(temp_c);

		if (dataset.size() == dataset_size) {
			break;
		}
	}

	fclose(input);

	return 0;
}


double LSH_clustering(const int& max_vec, const int& dim, const int& table_size, vector<int>& clusters, int& k)
{
    r.reserve(max_vec*K);
    for(int i=0; i<max_vec*K; i++)
    {
        int rr = (rand() / (RAND_MAX + 1.0))*(MAX_R+1);
        r.push_back(rr);
    }

    dataset_vectors.resize(dataset.size()*L);

    for(int l=0; l<L; l++)
    {
        g.push_back(new grid(DELTA, dim, K));
        table.push_back(new LSH_table(table_size, K_VEC, max_vec*K));

        for(int i=0; i<dataset.size(); i++)
        {
            curve temp_g(i, dim);

            g[l]->hash_curve(*dataset[i], &temp_g);

            dataset_vectors[i+(l*dataset.size())] = new vector_curve;
            dataset_vectors[i+(l*dataset.size())]->curve_to_vector(temp_g);
            dataset_vectors[i+(l*dataset.size())]->src = dataset[i];

            table[l]->insert_entry(dataset_vectors[i+(l*dataset.size())]);
        }
    }

    clusters.resize(dataset.size());
    vector<curve*> centroids;

    vector<int> ids;

    int cluster_idx = 0;
    for (int i=0; i<table_size; i++) {
    	table[0]->get_data(i)->get_ids(ids);

    	if (ids.size() > 0) {
    		centroids.push_back(dataset[ids[0]]);

    	    for (int j=0; j<ids.size(); j++) {
    			clusters[ids[j]] = cluster_idx;
    		}
    		cluster_idx++;
    	}
    }

    mean_update(dataset, clusters, centroids);

    vector<double> s;
    silhouette(dataset, clusters, centroids, s);

    double silhouette_coefficient = 0.0;
    for(int i=0; i<s.size(); i++)
    	silhouette_coefficient += s[i];

    k = centroids.size();

    return silhouette_coefficient /= s.size();
}

void LSH_handler(const int& max_nodes, const int& dim, const int& table_size) {
	clock_t t = clock();

	vector<int> clusters;
	int k;

	double s = LSH_clustering(max_nodes, dim, table_size, clusters, k);

	t = clock() - t;

	double time = (((double)t)/CLOCKS_PER_SEC) * 1000;

	print_output_rmsd(clusters, k, s, time, true);
}

void parse_config(char* file)
{
    FILE *conf;
    if((conf = fopen(file, "r")) == NULL)
        exit(1);

    char buf[32];
    int val;

    while(fscanf(conf, "%s %d", buf, &val) > 0)
    {
        if(strcmp(buf, "number_of_clusters:") == 0)
            k = val;
        else if(strcmp(buf, "number_of_grid_curves:") == 0)
            K = val;
        else if(strcmp(buf, "number_of_hash_tables:") == 0)
            L = val;
    }

    fclose(conf);
}

int main(int argc, char *argv[]) {
	if(argc != 9) {
        printf("Usage: $./segments -i <input_csv_unsegmented> -o <output_csv_segmented> -c <clustering_config_file> -s <size_of_dataset>\n");
        return 1;
    }

    srand(time(NULL));

    char *input_file, *conf_file, *segments_file;
    int max_nodes;

    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i], "-i") == 0)
            input_file = argv[i+1];
        else if (strcmp(argv[i], "-c") == 0)
            conf_file = argv[i+1];
        else if (strcmp(argv[i], "-s") == 0)
        	dataset_size = atoi(argv[i+1]);
        else if (strcmp(argv[i], "-o") == 0)
        	segments_file = argv[i+1];
    }

	srand(time(NULL));

	find_intersections(input_file);

	clock_t t = clock();
	if(parse(input_file, segments_file) != 0)
		return 1;
	t = clock() - t;

	printf("unsegmented -> segmented file parse time %f\n", ((double)t)/CLOCKS_PER_SEC);

	t = clock();
	if (parse_segments(segments_file, max_nodes) != 0) {
		return 1;
	}
	t = clock() - t;

	printf("segmented file parse time %f\n", ((double)t)/CLOCKS_PER_SEC);

    t = clock();
    parse_config(conf_file);
    t = clock() - t;

    printf("config file parsing time: %f\n", ((double)t)/CLOCKS_PER_SEC);

    L = 1;

    /*global initialization*/
    cached_distances.resize(dataset.size());
    for(int i=0; i<cached_distances.size(); i++)
        cached_distances[i].resize(dataset.size(), -1.0);
    /*global initialization*/

    metric = "classic_frechet";
    char *output_file = "output/lsh_ways_frechet.dat";
    if((output = fopen(output_file, "w")) == NULL)
    	return -1;

    t = clock();
    LSH_handler(max_nodes, dataset[0]->dim, dataset.size());
    t = clock() - t;

    printf("LSH clustering complete in %f seconds\n", ((double)t)/CLOCKS_PER_SEC);

    fclose(output);

    for(int i=0; i<cached_distances.size(); i++) {
    	for(int j=0; j<cached_distances[i].size(); j++) {
    		cached_distances[i][j] = -1.0;
    	}
    }

    translate_dataset_zero(dataset);

    metric = "DFD";
    output_file = "output/kmeans_ways_frechet.dat";
    if((output = fopen(output_file, "w")) == NULL)
    	return -1;

    t = clock();
    cluster_with_eval(I1A1U2_clustering, dataset, k, true);
    t = clock() - t;

    printf("kmeans clustering complete in %f seconds\n", ((double)t)/CLOCKS_PER_SEC);

    fclose(output);

    for(int i=0; i<dataset.size(); i++)
    	delete dataset[i];
    for(int i=0; i<dataset_vectors.size(); i++)
    	delete dataset_vectors[i];

    for(int i=0; i<L; i++)
    {
    	delete g[i];
    	delete table[i];
    }

	return 0;
}