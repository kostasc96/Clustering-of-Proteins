#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstring>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>

#include "curve.h"
#include "distance_metrics.h"
#include "clustering.h"

using namespace std;

#define inf 9999.0

double euclid_norm(const coord& a, const coord& b)
{
    double dist = 0.0;
    for(int i=0; i<a.dim; i++)
        dist += (a.c[i]-b.c[i])*(a.c[i]-b.c[i]);

    return sqrt(dist);
}

double discrete_frechet_distance(const curve& p, const curve& q, curve& mean, const bool& calc_mean) //DFD
{
    double** c = (double**)malloc(sizeof(double*)*p.length);

    int i, j;
    for(i=0; i<p.length; i++)
        c[i] = (double*)malloc(sizeof(double)*q.length);

    c[0][0] = euclid_norm(p.data[0], q.data[0]);
    for(j=1; j<q.length; j++)
    {
        double new_dist = euclid_norm(p.data[0], q.data[j]);
        c[0][j] = c[0][j-1] >= new_dist ? c[0][j-1] : new_dist;
    }
    for(i=1; i<p.length; i++)
    {
        double new_dist = euclid_norm(p.data[i], q.data[0]);
        c[i][0] = c[i-1][0] >= new_dist ? c[i-1][0] : new_dist;
    }

    for(i=1; i<p.length; i++)
    {
        for(j=1; j<q.length; j++)
        {
            double min_val = c[i-1][j];
            if(c[i-1][j-1] < min_val) min_val = c[i-1][j-1];
            if(c[i][j-1] < min_val) min_val = c[i][j-1];

            double new_dist = euclid_norm(p.data[i], q.data[j]);
            c[i][j] = min_val >= new_dist ? min_val : new_dist;
        }
    }

    double ret = c[p.length-1][q.length-1];

    if(calc_mean == true)
    {
        int p_idx = p.length-1, q_idx = q.length-1, min_val;

        coord temp_c;
        temp_c.c.resize(p.dim);
        temp_c.dim = p.dim;

        for(int i=0; i<p.dim; i++)
            temp_c.c[i] = (p.data[p_idx].c[i]+q.data[q_idx].c[i])/2;

        mean.data.reserve((p_idx+q_idx)/2);
        mean.data.insert(mean.data.begin(), temp_c);
        mean.length++;

        int max_len = (p.length > q.length) ? p.length : q.length;

        while((p_idx != 0 || q_idx != 0) && (mean.length < max_len))
        {
            if(p_idx == 0)
            {
                for(int i=0; i<p.dim; i++)
                    temp_c.c[i] = (p.data[0].c[i]+q.data[q_idx-1].c[i])/2;

                q_idx--;
            }
            else if(q_idx == 0)
            {
                for(int i=0; i<p.dim; i++)
                    temp_c.c[i] = (p.data[p_idx-1].c[i]+q.data[0].c[i])/2;

                p_idx--;
            }
            else
            {
                min_val = c[p_idx-1][q_idx-1];
                if(c[p_idx-1][q_idx] < min_val) min_val = c[p_idx-1][q_idx];
                if(c[p_idx][q_idx-1] < min_val) min_val = c[p_idx][q_idx-1];


                if(min_val == c[p_idx-1][q_idx-1])
                {
                    for(int i=0; i<p.dim; i++)
                        temp_c.c[i] = (p.data[p_idx-1].c[i]+q.data[q_idx-1].c[i])/2;

                    p_idx--;
                    q_idx--;
                }
                else if(min_val == c[p_idx-1][q_idx])
                {
                    for(int i=0; i<p.dim; i++)
                        temp_c.c[i] = (p.data[p_idx-1].c[i]+q.data[q_idx].c[i])/2;

                    p_idx--;
                }
                else
                {
                    for(int i=0; i<p.dim; i++)
                        temp_c.c[i] = (p.data[p_idx].c[i]+q.data[q_idx-1].c[i])/2;

                    q_idx--;
                }
            }

            mean.data.insert(mean.data.begin(), temp_c);
            mean.length++;
        }
    }

    for(i=0; i<p.length; i++)
        free(c[i]);
    free(c);

    return ret;
}

double dynamic_time_warping(const curve& p, const curve& q) //DTW
{
    double** c = (double**)malloc(sizeof(double*)*(p.length+1));

    int i, j;
    for(i=0; i<p.length+1; i++)
        c[i] = (double*)malloc(sizeof(double)*(q.length+1));

    c[0][0] = 0.0;
    for(j=1; j<q.length+1; j++)
        c[0][j] = inf;
    for(i=1; i<p.length+1; i++)
        c[i][0] = inf;

    for(i=1; i<p.length+1; i++)
    {
        for(j=1; j<q.length+1; j++)
        {
            double min_val = c[i-1][j];
            if(c[i-1][j-1] < min_val) min_val = c[i-1][j-1];
            if(c[i][j-1] < min_val) min_val = c[i][j-1];
            c[i][j] = euclid_norm(p.data[i-1], q.data[j-1]) + min_val;
        }
    }

    double ret = c[p.length][q.length];

    for(i=0; i<p.length+1; i++)
        free(c[i]);
    free(c);

    return ret;
}

using Eigen::MatrixXd;
using Eigen::JacobiSVD;

double c_RMSD(const curve& x, const curve& y)
{
    if (x.dim != y.dim) {
        return -1.0;
    }

    int dim = x.dim, length;

    curve x_t(-1, dim), y_t(-1, dim);

    if (x.length == y.length) {
        length = x.length;
    }
    else {
        if (x.length < y.length) {
            length = x.length;
            for (int i=0; i<x.length; i++) {
                y_t.add_coord(y.data[i]);
            }
        } else {
            length = y.length;
            for (int i=0; i<y.length; i++) {
                x_t.add_coord(x.data[i]);
            }
        }
    }

    MatrixXd X(length, dim);
    MatrixXd Y(length, dim);

    for (int i=0; i<length; i++) {
        for (int j=0; j<dim; j++) {
            X(i,j) = (x.length == length) ? x.data[i].c[j] : x_t.data[i].c[j];
            Y(i,j) = (y.length == length) ? y.data[i].c[j] : y_t.data[i].c[j];
        }
    }

    JacobiSVD<MatrixXd> svd(X.transpose()*Y, Eigen::ComputeFullV | Eigen::ComputeFullU);

    if(svd.singularValues()(dim-1) < 0)
        return -1.0;

    MatrixXd Q = svd.matrixU() * svd.matrixV().transpose();

    if(Q.determinant() < 0)
    {
        MatrixXd U = svd.matrixU();

        for (int i=0; i<dim; i++) {
            U(i, dim-1) = -U(i, dim-1);
        }

        Q = U * svd.matrixV().transpose();
    }

    if(strcmp(metric, "crmsd") == 0) {
        MatrixXd M = X * Q - Y;
        return sqrt((M.transpose()*M).trace())/sqrt(x.length);       
    }
    else {
        MatrixXd XQ = X * Q;

        curve xq(-1, dim);

        coord c;
        c.dim = dim;
        c.c.resize(dim);

        for (int i=0; i<XQ.rows(); i++) {
            for (int j=0; j<dim; j++) {
                c.c[j] = XQ(i, j);
            }

            xq.add_coord(c);
        }

        if(strcmp(metric, "DFD") == 0) {
            return discrete_frechet_distance(xq, ((y.length == length) ? y : y_t), *(curve*)NULL, false);
        } else if (strcmp(metric, "DTW") == 0) {
            return dynamic_time_warping(xq, ((y.length == length) ? y : y_t));
        }
    }
}

double calc_dist_update_cache(const curve& p, const int& p_idx, const curve& q, const int& q_idx)
{
    double dist;

    if((p_idx != -1 && q_idx != -1) && cached_distances[p_idx][q_idx] != -1.0)
        dist = cached_distances[p_idx][q_idx];
    else
    {
        if (!strcmp(metric, "classic_frechet")) {
            dist = discrete_frechet_distance(p, q, *(curve*)NULL, false);
        } else if (!strcmp(metric, "classic_DTW")) {
            dist = dynamic_time_warping(p, q);
        } else {
            dist = c_RMSD(p, q);            
        }

        if(p_idx != -1 && q_idx != -1)
            cached_distances[p_idx][q_idx] = cached_distances[q_idx][p_idx] = dist;
    }

    return dist;
}