#ifndef _CURVE_H_
#define _CURVE_H_

#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

extern FILE* files;

struct coord
{
    vector<double> c;
    int dim;

    bool cmp(const coord&);
    void print();
    void print_out();
    void subtract(const coord&);
    void add(const coord&);
    void divide(int);
};

class curve
{
    public:
        curve(int, int);
        ~curve();

        void add_coord(const coord&);
        bool cmp(const curve&);

        void print();

        void translate_zero();

        int id;

        int dim;
        int length;
        vector<coord> data;
};

struct vector_curve
{
    void curve_to_vector(const curve&);
    bool cmp(const vector_curve&);
    bool if_empty();
    void print();

    vector<double> data;

    curve* src;
    int dim;
    int id;
};

double dot_product(const vector<double>&, const vector<double>&, const int&);

void translate_dataset_zero(vector<curve*>&);

#endif // _CURVE_H_
