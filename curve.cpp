#include <vector>
#include <cstdio>

#include "curve.h"
#include "clustering.h"

bool coord::cmp(const coord& crd)
{
    for(int i=0; i<this->dim; i++)
    {
        if(this->c[i] != crd.c[i])
            return false;
    }
    return true;
}

void coord::print()
{
    fprintf(files, "(");
    for(int i=0; i<dim; i++)
    {
        fprintf(files, "%f", c[i]);
        if(i < dim-1)
            fprintf(files, ",");
    }
    fprintf(files, ")");
}

void coord::print_out()
{
    fprintf(output, "(");
    for(int i=0; i<dim; i++)
    {
        fprintf(output, "%f", c[i]);
        if(i < dim-1)
            fprintf(output, ",");
    }
    fprintf(output, ")");
}

void coord::subtract(const coord& c)
{
    if(this->dim != c.dim)
        return;

    for(int i=0; i<this->dim; i++)
        this->c[i] -= c.c[i];
}

void coord::add(const coord& c)
{
    if(this->dim != c.dim)
        return;

    for(int i=0; i<this->dim; i++)
        this->c[i] += c.c[i];
}

void coord::divide(int m)
{
    for(int i=0; i<this->dim; i++)
        this->c[i] /= (double)m;
}

curve::curve(int id, int dim):id(id), dim(dim)
{
    this->length = 0;
}

curve::~curve()
{

}

void curve::add_coord(const coord& c)
{
    data.push_back(c);
    length++;
}

bool curve::cmp(const curve& p)
{
    if(this->dim != p.dim || this->length != p.length)
        return false;

    for(int i=0; i<this->length; i++)
        if(data[i].cmp(p.data[i]) == false)
            return false;

    return true;
}

void curve::print()
{
    for(int i=0; i<length; i++)
    {
        data[i].print();
        if(i < length-1)
            fprintf(files, ",");
    }
    fprintf(files, "\n");
}

void curve::translate_zero()
{
    coord c;

    c.dim = this->dim;
    c.c.resize(this->dim);

    for (int i=0; i<this->dim; i++) {
        c.c[i] = 0.0;
    }

    for(int i=0; i<this->length; i++)
        c.add(this->data[i]);

    c.divide(this->length);

    for(int i=0; i<this->length; i++)
        this->data[i].subtract(c);
}

void vector_curve::curve_to_vector(const curve& c)
{
    if(c.length == 0)
        return;

    this->id = c.id;
    this->src = NULL;
    this->dim = 0;

    for(int i=0; i<c.length; i++)
    {
        for(int j=0; j<c.dim; j++)
        {
            data.push_back(c.data[i].c[j]);
            dim++;
        }
    }
}

bool vector_curve::cmp(const vector_curve& p)
{
    if(this->dim != p.dim)
        return false;

    for(int i=0; i<this->dim; i++)
    {
        if(this->data[i] != p.data[i])
            return false;
    }

    return true;
}

bool vector_curve::if_empty()
{
    return data.size() == 0;
}

void vector_curve::print()
{
    fprintf(files, "<");
    for(int i=0; i<dim; i++)
    {
        fprintf(files, "%f", data[i]);
        if(i < dim-1)
            fprintf(files, ",");
    }
    fprintf(files, ">\n");
}

double dot_product(const vector<double>& p, const vector<double>& q, const int& dim)
{
    double product = 0.0;

    for(int i=0; i<dim; i++)
        product += p[i]*q[i];

    return product;
}

void translate_dataset_zero(vector<curve*>& dataset)
{
    for(int i=0; i<dataset.size(); i++)
        dataset[i]->translate_zero();
}