#include <math.h>
#include <vector>
#include <algorithm>
#include "../utils.hpp"

namespace dF
{
    // recursive computation of table values for discrete Frechet distance
    double compute_c(vector<vector<double>> &c, int i, int j, const Item &p, const Item &q)
    {
        if (c[i][j] > -1)
        {
            return c[i][j];
        }
        else if ((i == 0) && (j > 0))
        {
            c[0][j] = max(compute_c(c, 0, j - 1, p, q), L2(0, p.xij[0], j, q.xij[j]));
        }
        else if ((i > 0) && (j == 0))
        {
            c[i][0] = max(compute_c(c, i - 1, 0, p, q), L2(i, p.xij[i], 0, q.xij[0]));
        }
        else if ((i > 0) && (j > 0))
        {
            double c1 = compute_c(c, i - 1, j, p, q);
            double c2 = compute_c(c, i - 1, j - 1, p, q);
            double c3 = compute_c(c, i, j - 1, p, q);
            double temp_min = min(c1, c2);
            double final_min = min(temp_min, c3);
            c[i][j] = max(final_min, L2(i, p.xij[i], j, q.xij[j]));
        }
        return c[i][j];
    }

    double discrete_frechet(const Item &p, const Item &q)
    {
        int pl = p.xij.size();
        int ql = q.xij.size();
        if (pl != ql)
        {
            std::cout << "Cannot compute discrete Frechet distance for curves with different dimensions." << std::endl;
            return -1;
        }
        else
        {
            // initialize 2d table for dynamic programming with value -1
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
            c[0][0] = L2(0, p.xij[0], 0, q.xij[0]);
            // begin computing
            return compute_c(c, pl - 1, ql - 1, p, q);
        }
    }

}
