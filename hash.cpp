#include <cstdlib>
#include <climits>
#include <cstdio>
#include <string>
#include <limits>
#include <cmath>

#include "distance_metrics.h"
#include "hash.h"

#define M 4294967291LL
#define W 4

using namespace std;

/*########## Normal Distribution #########*/

double randn(double mu, double sigma)
{
    double U1, U2, Z, mult;
    static double X1, X2;
    static int call = 0;

    if(call == 1)
    {
        call = !call;
        return (mu + sigma * (double)X2);
    }
    do
    {
        U1 = -1 + ((double)rand()/ RAND_MAX)*2;
        U2 = -1 + ((double)rand()/ RAND_MAX)*2;
        Z = (U1*U1) + (U2*U2);
    }
    while(Z >= 1 || Z == 0);

    mult = sqrt((-2 * log(Z))/ Z);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double)X1);
}

/*########## Normal Distribution #########*/

void gen_rand_vector(const int& dim, vector<double>& v)
{
    v.reserve(dim);

    for(int i=0; i<dim; i++)
    {
        double r = randn(0, 1);
        v.push_back(r);
    }
}

long long mod(long long k, const long long& n) {
    return ((k %= n) < 0) ? k+n : k;
}

/*########## curve list stuff ###########*/

curve_list::curve_list()
{
    head = NULL;
    tail = NULL;
    list_size = 0;
}

curve_list::~curve_list()
{
    node *current_node = head;
    node *temp_node = NULL;
	while(current_node != NULL)
	{
		temp_node = current_node;
		current_node = current_node->next;
		delete temp_node;
	}
}

void curve_list::insert_node(vector_curve* p)
{
    node *new_node = new node;
    new_node->data = p;
    new_node->next = NULL;

    if(head == NULL)
    {
        head = new_node;
        tail = new_node;
        list_size++;
    }
    else
    {
        tail->next = new_node;
        tail = tail->next;
        list_size++;
    }
}

bool curve_list::find_nn_r(vector_curve* p, const double& R, const int& c_id, vector<int>& ids, vector<int>& flagged, vector<int>& conflicts, double& of)
{
    node* current_node = head;

    if(list_size == 0)
        return false;

    double temp_dist;
    bool found = false;

    for(int i=0; i<list_size; i++)
    {
        if(flagged[current_node->data->id] != 2 && current_node->data->cmp(*p))
        {
            temp_dist = calc_dist_update_cache(*(p->src), p->id, *(current_node->data->src), current_node->data->id);

            if(temp_dist <= R)
            {
                if(flagged[current_node->data->id] == 0)
                {
                    ids[current_node->data->id] = c_id;
                    of += temp_dist;
                }
                else if(flagged[current_node->data->id] == 1)
                    conflicts.push_back(current_node->data->id);

                flagged[current_node->data->id]++;
                found = true;
            }
        }

        current_node = current_node->next;
    }

    if(found == true)
        return true;
    else
    {
        current_node = head;

        for(int i=0; i<list_size; i++)
        {
            if(flagged[current_node->data->id] != 2)
            {
                temp_dist = calc_dist_update_cache(*(p->src), p->id, *(current_node->data->src), current_node->data->id);

                if(temp_dist <= R)
                {
                    if(flagged[current_node->data->id] == 0)
                    {
                        ids[current_node->data->id] = c_id;
                        of += temp_dist;
                    }
                    else if(flagged[current_node->data->id] == 1)
                        conflicts.push_back(current_node->data->id);

                    flagged[current_node->data->id]++;
                    found = true;
                }
            }
 
            current_node = current_node->next;
        }

        if(found == true)
            return true;
    }

    return false;
}

void curve_list::print()
{
    node* current_node = head;
    for(int i=0; i<list_size; i++)
    {
        fprintf(files, "%d object: ", i);
        current_node->data->print();

        current_node = current_node->next;
    }
}

void curve_list::print_ids() {
    node* current_node = head;
    for (int i=0; i<list_size; i++) {
        printf("%d ", current_node->data->id);

        current_node = current_node->next;
    }
}

void curve_list::get_ids(vector<int>& ids) {
    ids.clear();

    node* current_node = head;
    for (int i=0; i<list_size; i++) {
        ids.push_back(current_node->data->id);

        current_node = current_node->next;
    }
}

/*########## curve list stuff ###########*/

/*########## hash table stuff ###########*/

classic_table::classic_table(int table_size):table_size(table_size)
{
    data = new curve_list[table_size];
}

classic_table::~classic_table()
{
    delete[] data;
}

LSH_table::LSH_table(int table_size, int k_vec, int max_dim):table_size(table_size), k_vec(k_vec)
{
    data = new curve_list[table_size];

    for(int i=0; i<k_vec; i++)
    {
        vector<double> temp;
        gen_rand_vector(max_dim, temp);
        v.push_back(temp);

        double rnd = (rand() / (RAND_MAX + 1.0))*W;
        rnd_t.push_back(rnd);
    }
}

LSH_table::~LSH_table()
{
    delete[] data;
}

void classic_table::insert_entry(vector_curve* p)
{
    int hash_index = h(p);
    data[hash_index].insert_node(p);
}

void LSH_table::insert_entry(vector_curve* p)
{
    int hash_index = h(p);
    data[hash_index].insert_node(p);
}

curve_list* LSH_table::get_data(int i) {
    return &data[i];
}

void classic_table::print()
{
    fprintf(files, "Printing hash table:\n");
    for(int i=0; i<table_size; i++)
    {
        fprintf(files, "Printing bucket %d:\n", i);
        data[i].print();
    }
}

void LSH_table::print()
{
    fprintf(files, "Printing hash table:\n");
    for(int i=0; i<table_size; i++)
    {
        fprintf(files, "Printing bucket %d:\n", i);
        data[i].print();
    }
}

/*################ hash functions #################*/
int classic_table::h(vector_curve* p)
{
    long long lc = 0;

    for(int i=0; i<p->dim; i++)
        lc = mod((lc + mod((mod(r[i],M) * mod((int)p->data[i],M)),M)),M);

    lc = mod(lc,M);

    return mod(lc,table_size);
}

int LSH_table::h(vector_curve* p)
{
    long long f = 0;

    for(int i=0; i<k_vec; i++)
    {
        int h = (int)((dot_product(p->data,v[i],p->dim) + rnd_t[i])/(double)W);

        f = mod((f + mod((mod(r[i],M) * mod(h,M)),M)),M);
    }

    f = mod(f,M);

    return mod(f,table_size);
}
/*################ hash functions #################*/

bool classic_table::find_nn_r(vector_curve* p, const double& R, const int& c_id, vector<int>& ids, vector<int>& flagged, vector<int>& conflicts, double& of)
{
    int index = h(p);
    return data[index].find_nn_r(p, R, c_id, ids, flagged, conflicts, of);
}

bool LSH_table::find_nn_r(vector_curve* p, const double& R, const int& c_id, vector<int>& ids, vector<int>& flagged, vector<int>& conflicts, double& of)
{
    int index = h(p);
    return data[index].find_nn_r(p, R, c_id, ids, flagged, conflicts, of);
}

/*########## hash table stuff ###########*/