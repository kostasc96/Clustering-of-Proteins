#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>
#include <cmath>

#include "curve.h"
#include "grid.h"
#include "hash.h"
#include "clustering.h"
#include "assignment.h"
#include "update.h"
#include "tree.h"

using namespace std;

#define DELTA 0.02
#define MAX_R 100

void parse_config(char*);

void parse_file(char*);
void create_hashtables(const int&, const int&, const int&);

int L;
int K;
int k;

vector< vector<double> > cached_distances;
vector<int> r;

vector<curve*> dataset;
vector<vector_curve*> dataset_vectors;

vector<grid*> g;
vector<LSH_table*> table;

clock_t t;

FILE *output;
FILE *files;

char *metric;

int main(int argc, char* argv[])
{
    if (argc != 5) {
        printf("Usage: $./proteins -i <bio_input_file> -c <clustering_config_file>\n");
        return 1;
    }

    srand(time(NULL));

    char *input_file, *conf_file, *output_file;
    int max_vec = 0, table_size;

    for(int i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "-i") == 0)
            input_file = argv[i+1];
        else if(strcmp(argv[i], "-c") == 0)
            conf_file = argv[i+1];
    }

    t = clock();
    parse_config(conf_file);
    t = clock() - t;

    printf("config file parsing time: %f\n", ((double)t)/CLOCKS_PER_SEC);

    /*input*/
    t = clock();
    parse_file(input_file);
    t = clock() - t;

    printf("input file parsing time: %f\n", ((double)t)/CLOCKS_PER_SEC);
    /*input*/

    /*global initialization*/
    cached_distances.resize(dataset.size());
    for(int i=0; i<cached_distances.size(); i++)
        cached_distances[i].resize(dataset.size(), -1.0);

    translate_dataset_zero(dataset);
    /*global initialization*/

    metric = "crmsd";
    output_file = "output/crmsd.dat";
    if((output = fopen(output_file, "w")) == NULL)
    	return -1;

    t = clock();
    cluster_with_eval(I1A1U2_clustering, dataset, k, false);
    t = clock() - t;

    printf("crmsd clustering complete in %f seconds\n", ((double)t)/CLOCKS_PER_SEC);

    fclose(output);

    for(int i=0; i<cached_distances.size(); i++) {
    	for(int j=0; j<cached_distances[i].size(); j++) {
    		cached_distances[i][j] = -1.0;
    	}
    }

    metric = "DFD";
    output_file = "output/frechet.dat";
    if((output = fopen(output_file, "w")) == NULL)
    	return -1;

    t = clock();
    cluster_with_eval(I1A1U2_clustering, dataset, k, false);
    t = clock() - t;

    printf("frechet clustering complete in %f seconds\n", ((double)t)/CLOCKS_PER_SEC);

    fclose(output);

    for(int i=0; i<dataset.size(); i++)
    	delete dataset[i];
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

void parse_file(char* file) //used to parse the input
{
    FILE *input;
    if((input = fopen(file, "r")) == NULL)
        exit(1);

    int numConform, N;

    fscanf(input, "%d\n%d\n", &numConform, &N);

    coord c;
    c.dim = 3;
    c.c.resize(3);

    for(int i=0; i<numConform; i++)
    {
    	curve* temp_c = new curve(dataset.size(), 3);
        temp_c->data.reserve(N);

        for(int j=0; j<N; j++)
        {
        	fscanf(input, "%lf %lf %lf\n", &(c.c[0]), &(c.c[1]), &(c.c[2]));

        	temp_c->add_coord(c);
        }

        dataset.push_back(temp_c);
    }

    fclose(input);
}

void create_hashtables(const int& max_vec, const int& dim, const int& table_size)
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
}