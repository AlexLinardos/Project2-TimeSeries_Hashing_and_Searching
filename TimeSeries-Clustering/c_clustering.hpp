#ifndef C_CLUSTERING_HPP
#define C_CLUSTERING_HPP

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include "../includes/curves.hpp"
#include "../includes/utils.hpp"
#include "../includes/mean_curve.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/disc_Frechet.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/discF_LSH.hpp"
#include "../TimeSeries-ANN/L2/LSH.hpp"
#include "../TimeSeries-ANN/L2/HC.hpp"

namespace curve_cluster
{
    class Clustering
    {
        int n_centers;  // number of centers to be initialized
        std::vector<curves::Curve2d> * dataset; // a dataset of curves

        vector<vector<curves::Point2d>> centers;
        vector<vector<curves::Curve2d>> clusters;
        std::vector<int> assignments_vec;     // shows the cluster to which each curve is assigned to
                                              // (example: if assignments[4]=2 then curve at index 4 of dataset is assigned to cluster at index 2)

        std::default_random_engine eng;       // an engine to help us pick first center randomnly
        std::uniform_int_distribution<int> uid;

        /* Calculates the Frechet distance of each point in dataset from its nearest centroid
            Parameters are an std::vector containing the current centroids, and std::vector in which the min distance from nearest centroid
            of each corresponding dataset item will be stored and a sum accumulator "dist_sum" of all those distances (needed for calculate_probs).*/
        void calculate_min_dists(vector<double> &d, double &dist_sum)
        {
            double max_d = 0.0; // for normalization
            // for each point in the dataset
            for (int i = 0; i < this->dataset->size(); ++i)
            {
                // calculate Frechet distance from nearest centoid
                double min_dist = dF::discrete_frechet_for_data(centers[0], (*dataset)[i].data);
                if (min_dist > max_d)
                {
                    max_d = min_dist;
                }
                for (int cntr = 1; cntr < centers.size(); ++cntr)
                {
                    double dist = dF::discrete_frechet_for_data(centers[cntr], (*dataset)[i].data);
                    if (dist > max_d)
                        max_d = dist;
                    if (dist < min_dist)
                        min_dist = dist;
                }
                d[i] = min_dist;
            }
            // normalize data
            for (int j = 0; j < d.size(); ++j)
            {
                d[j] = d[j] / max_d;
            }
            // calculate sum ΣD(i)^2
            for (int i = 0; i < d.size(); ++i)
            {
                dist_sum += pow(d[i], 2.0);
            }
        }

        /* Calculates probability of each point in dataset to be picked as next centroid.
            Parameters are an std::vector containing all the corresponding distances of each point in the dataset
            from its nearest centroid, a sum of all those distances and an std::vector in which the corresponding
            probability of each point will be stored */
        void calculate_probs(vector<double> d, double dist_sum, vector<double> &p)
        {
            for (int i = 0; i < d.size(); ++i)
            {
                p[i] = pow(d[i], 2.0) / dist_sum;
            }
        }

        // helper method that calculates starting radius for Reverse Assignment algorithms, which is min(dist between centers)/2
        double calculate_start_radius()
        {
            double min_dist = std::numeric_limits<double>::max();
            for (int i = 0; i < this->centers.size(); ++i)
            {
                for (int j = 0; j < this->centers.size(); ++j)
                {
                    if (i != j)
                    {
                        double dist = dF::discrete_frechet_for_data(this->centers[i], this->centers[j]);
                        if (dist < min_dist)
                        {
                            min_dist = dist;
                        }
                    }
                }
            }
            return min_dist / 2;
        }

        public:
        Clustering(int n_centers, std::vector<curves::Curve2d> * dataset) : n_centers(n_centers),
                                                                                dataset(dataset),
                                                                                eng(chrono::system_clock::now().time_since_epoch().count()),
                                                                                uid(0, dataset->size() - 1) {}

        void initialize_pp()
        {
            // pick first centroid at random
            int rcentroid_index = this->uid(eng);
            centers.push_back((*dataset)[rcentroid_index].data);

            // continue by picking more centroids until we reach the given number
            // centroids must be far away from each other
            while (centers.size() < this->n_centers)
            {
                double dist_sum = 0.0;
                std::vector<double> d(dataset->size(), 0.0); // distances of points to closest centroid
                std::vector<double> p(dataset->size(), 0.0); // probabilities of points to be chosen as next centroid

                this->calculate_min_dists(d, dist_sum); // calculate min distance of all points to the closest centroid for them
                // cout << "SUM (after normalization) = " << dist_sum << endl;

                this->calculate_probs(d, dist_sum, p); // calculate probability of each point to be chosen as next centroid
                double sum_prob = 0.0;
                for (int k = 0; k < p.size(); ++k)
                {
                    sum_prob += p[k];
                }

                std::discrete_distribution<> distribution{p.begin(), p.end()};
                int pick = distribution(eng);
                centers.push_back((*dataset)[pick].data);
            }
        }
        /* Assigns a nearest center to each point (part 1 of Lloyd's algorithm) using L2 distance.
         Stores both the current assignments and the current clusters in "assignments_vec" and "clusters" attributes respectively */
        void Lloyds_assignment()
        {
            int nearest_cntr;
            for (int i = 0; i < this->dataset->size(); ++i)
            {

                if ((*dataset)[i].marked) // if item is marked, it is assigned in a cluster during reverse_assignement
                    continue;

                double min_d = dF::discrete_frechet_for_data(this->centers[0], (*dataset)[i].data);
                nearest_cntr = 0;

                for (int c = 1; c < centers.size(); ++c)
                {
                    double next_d = dF::discrete_frechet_for_data(this->centers[c], (*dataset)[i].data);
                    if (next_d < min_d)
                    {
                        min_d = next_d;
                        nearest_cntr = c;
                    }
                }

                this->assignments_vec[i] = nearest_cntr;
                this->clusters[nearest_cntr].push_back((*dataset)[i]); // push it into a cluster based on assigned center
                (*dataset)[i].cluster = nearest_cntr;
            }
        }

        /* Assigns a nearest center to each point by reverse LSH range search and discrete Frechet distance.
         Expects a dFLSH::LSH object as parameter. */
        void Range_dfLSH_assignment(dFLSH::LSH &dflsh)
        {
            int balls_changed = 0;
            double radius = calculate_start_radius();
            int iter = 0;

            unordered_map<string, pair<curves::Curve2d *, int>> step_assignments; // map item id string to a pair that indicates a temporary item-cluster assignment
            do
            {
                balls_changed = 0;
                for (int c = 0; c < this->n_centers; ++c)
                {
                    // perform LSH Range Search
                    std::vector<std::pair<curves::Curve2d *, double>> r_search;
                    r_search = dflsh.RangeSearch(this->centers[c], radius, 0);

                    // if ball found new items
                    if (r_search.size() > 0)
                    {
                        balls_changed++;
                    }
                    for (int i = 0; i < r_search.size(); ++i)
                    {
                        curves::Curve2d *item = r_search[i].first;
                        if (item->claimed == false) // if item has not been claimed by a ball
                        {
                            item->claimed = true;

                            step_assignments[item->id] = make_pair(item, c); // temp assignment of item to cluster of index c
                        }
                        else // else must resolve conflict
                        {
                            // if we are here the item has already been assigned to another cluster
                            int assigned_cluster = step_assignments[item->id].second;

                            double dist_to_assigned = dF::discrete_frechet_for_data(this->centers[assigned_cluster], (*item).data);

                            double dist_to_curr = dF::discrete_frechet_for_data(this->centers[c], (*item).data);

                            if (dist_to_curr < dist_to_assigned)
                                step_assignments[item->id] = make_pair(item, c); // temp assignment of item to cluster of index c
                        }
                    }
                    r_search.clear();
                }
                radius *= 2; // double the search radius
                iter++;

                for (auto x : step_assignments)
                {
                    // get the closest cluster found
                    int closest_cluster = x.second.second;
                    // get the item
                    curves::Curve2d *item = x.second.first;

                    item->claimed = false;
                    // mark the item because it will be added in a cluster
                    item->marked = true;
                    // assign point to its closest cluster
                    item->cluster = closest_cluster;
                    clusters[closest_cluster].push_back(*item);
                }
                step_assignments.clear();

            } while (balls_changed >= this->n_centers * 0.2 || iter < 5);

            this->Lloyds_assignment();
        }

        // update step of clustering for curves by calculation of mean curve
        void update_curve_centers()
        {
            // we must calculate mean per cluster and make it the new center
            for (int i = 0; i < centers.size(); ++i)
            {
                centers[i] = mean_of_curves(clusters[i]);
            }
        }
    };
}
#endif