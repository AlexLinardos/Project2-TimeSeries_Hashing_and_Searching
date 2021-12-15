#ifndef DISC_FRECHET_HPP
#define DISC_FRECHET_HPP
#include <math.h>
#include <vector>
#include <algorithm>
#include "../utils.hpp"
#include "../curves.hpp"

namespace dF
{
    // recursive computation of table values for discrete Frechet distance (dynamic programming)
    double compute_c(vector<vector<double>> &c, int i, int j, const curves::Curve2d &p, const curves::Curve2d &q)
    {
        // if this value has already been computed return it immediately
        if (c[i][j] > -1)
        {
            return c[i][j];
        }
        else if ((i == 0) && (j > 0))
        {
            c[0][j] = std::max(compute_c(c, 0, j - 1, p, q), L2(p.data[0].x, p.data[0].y, q.data[j].x, q.data[j].y));
        }
        else if ((i > 0) && (j == 0))
        {
            c[i][0] = std::max(compute_c(c, i - 1, 0, p, q), L2(p.data[i].x, p.data[i].y, q.data[0].x, q.data[0].y));
        }
        else if ((i > 0) && (j > 0))
        {
            double c1 = compute_c(c, i - 1, j, p, q);
            double c2 = compute_c(c, i - 1, j - 1, p, q);
            double c3 = compute_c(c, i, j - 1, p, q);
            double temp_min = std::min(c1, c2);
            double final_min = std::min(temp_min, c3);
            c[i][j] = std::max(final_min, L2(p.data[i].x, p.data[i].y, q.data[j].x, q.data[j].y));
        }
        return c[i][j];
    }

    double discrete_frechet(const curves::Curve2d &p, const curves::Curve2d &q)
    {
        int pl = p.data.size();
        int ql = q.data.size();
        if (pl != ql)
        {
            std::cout << "Cannot compute discrete Frechet distance for curves with different dimensions." << std::endl;
            return -1;
        }
        else
        {
            // initialize 2d table for dynamic programming with value -1 (since distance cannot be negative)
            vector<vector<double>> c(pl);
            for (int i = 0; i < pl; i++)
            {
                c[i] = vector<double>(ql);
                for (int j = 0; j < ql; j++)
                {
                    c[i][j] = -1;
                }
            }
            // base case
            c[0][0] = L2(p.data[0].x, p.data[0].y, q.data[0].x, q.data[0].y);
            // begin computing
            return compute_c(c, pl - 1, ql - 1, p, q);
        }
    }

}
#endif