#ifndef UPDATE_HPP
#define UPDATE_HPP
#include <iostream>
#include <string>
#include <vector>
#include "../includes/utils.hpp"
#include "../includes/mean_curve.hpp"

namespace updt
{
    // update step of clustering for vectors by calculation of mean vector
    void update_vector_centers(std::vector<Item> &dataset, std::vector<Item> &centers, std::vector<vector<Item>> clusters)
    {
        int v_dimension = dataset[0].xij.size();
        // we must calculate mean per cluster and make it the new center
        for (int i = 0; i < centers.size(); ++i)
        {
            vector<double> mean(v_dimension, 0.0);
            int T = clusters[i].size(); // number of objects in cluster

            // calculate mean
            for (int j = 0; j < T; ++j)
            {
                mean = vector_mean(mean, clusters[i][j].xij, v_dimension, T); // using vector_addition from utilities.hpp
            }

            centers[i].xij = mean;
        }
    }

    // update step of clustering for curves by calculation of mean curve
    void update_curve_centers(std::vector<curves::Curve2d> &dataset, std::vector<curves::Curve2d> &centers, std::vector<vector<curves::Curve2d>> clusters)
    {
        // we must calculate mean per cluster and make it the new center
        for (int i = 0; i < centers.size(); ++i)
        {
            // CODE HERE
        }
    }
}

#endif