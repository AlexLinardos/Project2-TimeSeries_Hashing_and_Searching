#ifndef DISC_FRECHET_HPP
#define DISC_FRECHET_HPP
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../includes/utils.hpp"
#include "../../includes/curves.hpp"

namespace dF
{
    // recursive computation of table values for discrete Frechet distance (dynamic programming)
    double compute_c(double **c, int i, int j, const std::vector<curves::Point2d> &p, const std::vector<curves::Point2d> &q)
    {
        // if this value has already been computed return it immediately
        if (c[i][j] > -1)
        {
            return c[i][j];
        }
        else if ((i == 0) && (j > 0))
        {
            c[0][j] = std::max(compute_c(c, 0, j - 1, p, q), L2(p[0].x, p[0].y, q[j].x, q[j].y));
        }
        else if ((i > 0) && (j == 0))
        {
            c[i][0] = std::max(compute_c(c, i - 1, 0, p, q), L2(p[i].x, p[i].y, q[0].x, q[0].y));
        }
        else if ((i > 0) && (j > 0))
        {
            double c1 = compute_c(c, i - 1, j, p, q);
            double c2 = compute_c(c, i - 1, j - 1, p, q);
            double c3 = compute_c(c, i, j - 1, p, q);
            double temp_min = std::min(c1, c2);
            double final_min = std::min(temp_min, c3);
            c[i][j] = std::max(final_min, L2(p[i].x, p[i].y, q[j].x, q[j].y));
        }
        return c[i][j];
    }

    // calculates discrete Frechet distance between two curves (returns table of dynamic programming approach)
    double **discrete_frechet_table(const curves::Curve2d &p, const curves::Curve2d &q)
    {
        int pl = p.data.size();
        int ql = q.data.size();

        double **c = new double *[pl];
        try
        {
            for (int i = 0; i < pl; i++)
            {
                c[i] = new double[ql];
            }

            // initialize 2d table for dynamic programming with value -1 (since distance cannot be negative)
            for (int i = 0; i < pl; i++)
            {
                for (int j = 0; j < ql; j++)
                {
                    c[i][j] = -1;
                }
            }
        }
        catch (...)
        {
            std::cout << "INITIALIZATION" << endl;
            std::cout << "pl: " << pl << " ql: " << ql << endl;
        }

        // base case
        c[0][0] = L2(p.data[0].x, p.data[0].y, q.data[0].x, q.data[0].y);

        // compute
        compute_c(c, pl - 1, ql - 1, p.data, q.data);

        return c;
    }

    // calculates discrete Frechet distance between two curves (returns only final value)
    double discrete_frechet(const curves::Curve2d &p, const curves::Curve2d &q)
    {
        int pl = p.data.size();
        int ql = q.data.size();

        double **c = new double *[pl];

        for (int i = 0; i < pl; i++)
        {
            c[i] = new double[ql];
        }

        // initialize 2d table for dynamic programming with value -1 (since distance cannot be negative)
        for (int i = 0; i < pl; i++)
        {
            for (int j = 0; j < ql; j++)
            {
                c[i][j] = -1;
            }
        }

        // base case
        c[0][0] = L2(p.data[0].x, p.data[0].y, q.data[0].x, q.data[0].y);

        // compute
        compute_c(c, pl - 1, ql - 1, p.data, q.data);

        double result = c[pl - 1][ql - 1];

        for (int i = 0; i < pl; i++)
        {
            delete[] c[i];
        }
        delete[] c;

        return result;
    }

    // calculates discrete Frechet distance between only the data attribute of two curves (returns only final value)
    double discrete_frechet_for_data(std::vector<curves::Point2d> &p, std::vector<curves::Point2d> &q)
    {
        int pl = p.size();
        int ql = q.size();

        double **c = new double *[pl];

        for (int i = 0; i < pl; i++)
        {
            c[i] = new double[ql];
        }

        // initialize 2d table for dynamic programming with value -1 (since distance cannot be negative)
        for (int i = 0; i < pl; i++)
        {
            for (int j = 0; j < ql; j++)
            {
                c[i][j] = -1;
            }
        }

        // base case
        c[0][0] = L2(p[0].x, p[0].y, q[0].x, q[0].y);

        // compute
        compute_c(c, pl - 1, ql - 1, p, q);

        double result = c[pl - 1][ql - 1];

        for (int i = 0; i < pl; i++)
        {
            delete[] c[i];
        }
        delete[] c;

        return result;
    }

    // searches for the exact nearest neighbour of the query curve using brute force approach
    std::pair<curves::Curve2d *, double> search_exactNN(curves::Curve2d &query, vector<curves::Curve2d> &dataset)
    {
        int size = query.data.size(); // length of curve

        // we will store current nearest neighbour in curr_NN along with its distance from query
        vector<double> dummy_vec;
        for (int i = 0; i < size; i++)
            dummy_vec.push_back(0.0);
        curves::Curve2d null_curve = curves::Curve2d("null", dummy_vec, dummy_vec);
        std::pair<curves::Curve2d *, double> curr_NN;
        curr_NN.first = &null_curve;
        curr_NN.second = std::numeric_limits<double>::max();

        // for each curve in the dataset
        for (int i = 0; i < dataset.size(); i++)
        {
            // calculate discrete Frechet distance to it from given query
            double dfd = dF::discrete_frechet(query, dataset[i]);
            // if nearer curve is found
            if (dfd < curr_NN.second)
            {
                // replace curr_NN
                curr_NN.first = &(dataset[i]);
                curr_NN.second = dfd;
            }
        }
        return curr_NN;
    }
}
#endif