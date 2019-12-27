#ifndef _HASH_H
#define _HASH_H

#include <vector>

#include "curve.h"

#define K_VEC 3

extern int L;

extern vector<int> r;

int classic_h(const vector_curve&, int&);
int LSH_h(const vector_curve&, const vector< vector<double> >&, const int&, const int&);

class curve_list
{
    public:
        curve_list();
        ~curve_list();

        void insert_node(vector_curve*);
        bool find_nn_r(vector_curve*, const double&, const int&, vector<int>&, vector<int>&, vector<int>&, double&);

        void print();
        void print_ids();

        void get_ids(vector<int>&);

    private:
        struct node
        {
            vector_curve* data;
            node* next;
        };

        node *head;
        node *tail;
        int list_size;
};

class classic_table
{
    public:
        classic_table(int);
        ~classic_table();

        int h(vector_curve*);
        bool find_nn_r(vector_curve*, const double&, const int&, vector<int>&, vector<int>&, vector<int>&, double&);
        void insert_entry(vector_curve*);

        void print();

    private:
        curve_list* data;
        int table_size;
};

class LSH_table
{
    public:
        LSH_table(int, int, int);
        ~LSH_table();

        int h(vector_curve*);
        bool find_nn_r(vector_curve*, const double&, const int&, vector<int>&, vector<int>&, vector<int>&, double&);
        void insert_entry(vector_curve*);

        curve_list* get_data(int);

        void print();

    private:
        curve_list* data;
        int table_size;

        int k_vec;
        vector< vector<double> > v;
        vector<double> rnd_t;
};

#endif // _HASH_H