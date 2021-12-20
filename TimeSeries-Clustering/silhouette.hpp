#ifndef SILHOUETTE_HPP
#define SILHOUETTE_HPP
#include <iostream>
#include <vector>
#include <string>
#include "../includes/utils.hpp"
#include "../includes/curves.hpp"

namespace sil
{
    // Silhouette of object at index i
    double silhouette(Item *item)
    {
        vector<Item> *cluster = &clusters[item->cluster]; // find out at which cluster this item is assigned to

        //  calculate a(i) = average distance of i to objects in same cluster
        vector<double> distances;

        for (int j = 0; j < (*cluster).size(); ++j)
        {
            double dist = EuclideanDistance(item, &((*cluster)[j]), dimension);
            distances.push_back(dist);
        }

        double a = 0.0;
        for (int j = 0; j < distances.size(); ++j)
        {
            a += distances[j] / (double)distances.size();
        }

        // calculate b(i) = average distance of i to objects in next best (neighbor) cluster, i.e. cluster of 2nd closest centroid
        // first find the index of the next best cluster
        int best;
        double best_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < this->centers.size(); ++j)
        {
            if (j != item->cluster)
            {
                double dist = EuclideanDistance(item, &this->centers[j], dimension);
                if (dist < best_dist)
                {
                    best_dist = dist;
                    best = j;
                }
            }
        }

        // now perform the calculations
        distances.clear();
        cluster = &(this->clusters[best]);
        for (int j = 0; j < (*cluster).size(); ++j)
        {
            double dist = EuclideanDistance(item, &((*cluster)[j]), dimension);
            distances.push_back(dist);
        }
        double b = 0.0;
        for (int j = 0; j < distances.size(); ++j)
        {
            b += distances[j] / (double)distances.size();
        }
        // find max between a(i) and b(i)
        double max = a;
        if (b > a)
        {
            max = b;
        }
        // finally calculate and return silhouette
        return (b - a) / max;
    }

    // metric to evaluate specific cluster (c_index = the index of said cluster)
    double eval_specific_cluster(int c_index)
    {
        // cout << "Evaluating cluster: " << c_index << endl;
        double average = 0.0;
        double s = 0.0;
        int n = this->clusters[c_index].size();
        for (int i = 0; i < n; i++)
        {
            s = silhouette(&clusters[c_index][i]);
            average += s / (double)n;
        }
        return average;
    }

    // metric to evaluate overall clustering
    double eval_clustering()
    {
        double average = 0.0;
        double s = 0.0;
        int n = this->dataset.size();
        for (int i = 0; i < n; ++i)
        {
            s = silhouette(&dataset[i]);
            average += s / (double)n;
        }
        return s;
    }
}
#endif