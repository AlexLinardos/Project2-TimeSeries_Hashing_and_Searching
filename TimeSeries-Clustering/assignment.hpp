#ifndef ASSIGNMENT_HPP
#define ASSIGNMENT_HPP
#include <iostream>
#include <vector>
#include <string>
#include "../includes/curves.hpp"
#include "../includes/utils.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/disc_Frechet.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/discF_LSH.hpp"
#include "../TimeSeries-ANN/L2/LSH.hpp"
#include "../TimeSeries-ANN/L2/HC.hpp"

namespace assign
{
    // helps us implement the "assignment" step of the clustering procedure for vectors with various algorithms (uses L2 distance)
    class VectorAssignor
    {
    private:
        int n_centers;                    // number of centers/clusters
        std::vector<Item> dataset;        // a dataset of vectors
        int dimensions;                   // dimensions of vectors in the dataset
        std::vector<int> assignments_vec; // shows the cluster to which each vector is assigned to
                                          // (example: if assignments[4]=2 then item at index 4 of dataset is assigned to cluster at index 2)

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
                        double dist = EuclideanDistance(&this->centers[i], &this->centers[j], this->dimensions);
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
        std::vector<Item> centers;          // cluster centers (constructor receives them initialized)
        std::vector<vector<Item>> clusters; // the clusters of items
        VectorAssignor(int n_centers, std::vector<Item> centers, std::vector<Item> dataset, int dimensions) : n_centers(n_centers),
                                                                                                              centers(centers),
                                                                                                              dataset(dataset),
                                                                                                              dimensions(dimensions),
                                                                                                              clusters(n_centers),
                                                                                                              assignments_vec(dataset.size()) {}

        /* Assigns a nearest center to each point (part 1 of Lloyd's algorithm) using L2 distance.
         Stores both the current assignments and the current clusters in "assignments_vec" and "clusters" attributes respectively */
        void Lloyds()
        {
            int nearest_cntr;
            for (int i = 0; i < this->dataset.size(); ++i)
            {
                if (dataset[i].marked) // if item is marked, it is assigned in a cluster during reverse_assignement
                    continue;

                double min_d = EuclideanDistance(&(centers[0]), &(this->dataset[i]), this->dimensions);
                nearest_cntr = 0;
                for (int c = 1; c < centers.size(); ++c)
                {
                    double next_d = EuclideanDistance(&centers[c], &this->dataset[i], this->dimensions);
                    if (next_d < min_d)
                    {
                        min_d = next_d;
                        nearest_cntr = c;
                    }
                }
                this->assignments_vec[i] = nearest_cntr;
                this->clusters[nearest_cntr].push_back(this->dataset[i]); // push it into a cluster based on assigned center
                dataset[i].cluster = nearest_cntr;
            }
        }

        /* Assigns a nearest center to each point by reverse LSH range search using L2 distance.
          Expects an LSH object as parameter. */
        void Range_LSH(LSH &lsh)
        {
            int balls_changed = 0;
            double radius = calculate_start_radius();
            int iter = 0;

            unordered_map<string, pair<Item *, int>> step_assignments; // map item id string to a pair that indicates a temporary item-cluster assignment
            do
            {
                balls_changed = 0;
                for (int c = 0; c < this->n_centers; ++c)
                {
                    // perform LSH Range Search
                    std::vector<std::pair<double, Item *>> r_search;
                    r_search = lsh.RangeSearch(&this->centers[c], radius, 0);

                    // if ball found new items
                    if (r_search.size() > 0)
                    {
                        balls_changed++;
                    }
                    for (int i = 0; i < r_search.size(); ++i)
                    {
                        Item *item = r_search[i].second;
                        if (item->claimed == false) // if item has not been claimed by a ball
                        {
                            item->claimed = true;

                            step_assignments[item->id] = make_pair(item, c); // temp assignment of item to cluster of index c
                        }
                        else // else must resolve conflict
                        {
                            // if we are here the item has already been assigned to another cluster
                            int assigned_cluster = step_assignments[item->id].second;

                            double dist_to_assigned = EuclideanDistance(&centers[assigned_cluster], item, this->dimensions);
                            double dist_to_curr = EuclideanDistance(&centers[c], item, this->dimensions);

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
                    Item *item = x.second.first;

                    item->claimed = false;
                    // mark the item because it will be added in a cluster
                    item->marked = true;
                    // assign point to its closest cluster
                    item->cluster = closest_cluster;
                    clusters[closest_cluster].push_back(*item);
                }
                step_assignments.clear();

            } while (balls_changed >= this->n_centers * 0.2 || iter < 5);

            this->Lloyds();
        }

        /* Assigns a nearest center to each point by reverse Hypercube range search using L2 distance.
          Expects a Hypercube object as parameter. */
        void Range_HC(Hypercube &hc)
        {
            int balls_changed = 0;
            double radius = calculate_start_radius();
            int iter = 0;

            unordered_map<string, pair<Item *, int>> step_assignments; // map item id string to a pair that indicates a temporary item-cluster assignment
            do
            {
                balls_changed = 0;
                for (int c = 0; c < this->n_centers; ++c)
                {
                    // perform Hypercube Range Search
                    std::vector<std::pair<double, Item *>> r_search;
                    r_search = hc.RangeSearch(&this->centers[c], radius);

                    // if ball found new items
                    if (r_search.size() > 0)
                    {
                        balls_changed++;
                    }
                    for (int i = 0; i < r_search.size(); ++i)
                    {
                        Item *item = r_search[i].second;
                        if (item->claimed == false) // if item has not been claimed by a ball
                        {
                            item->claimed = true;

                            step_assignments[item->id] = make_pair(item, c); // temp assignment of item to cluster of index c
                        }
                        else // else must resolve conflict
                        {
                            // if we are here the item has already been assigned to another cluster
                            int assigned_cluster = step_assignments[item->id].second;

                            double dist_to_assigned = EuclideanDistance(&centers[assigned_cluster], item, this->dimensions);
                            double dist_to_curr = EuclideanDistance(&centers[c], item, this->dimensions);

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
                    Item *item = x.second.first;

                    item->claimed = false;
                    // mark the item because it will be added in a cluster
                    item->marked = true;
                    // assign point to its closest cluster
                    item->cluster = closest_cluster;
                    clusters[closest_cluster].push_back(*item);
                }
                step_assignments.clear();

            } while (balls_changed >= this->n_centers * 0.2 || iter < 5);

            this->Lloyds();
        }
    };

    // helps us implement the "assignment" step of the clustering procedure for curves with various algorithms (uses discrete Frechet distance)
    class CurveAssignor
    {
    private:
        int n_centers;                        // number of centers/clusters
        std::vector<curves::Curve2d> dataset; // a dataset of curves
        std::vector<int> assignments_vec;     // shows the cluster to which each curve is assigned to
                                              // (example: if assignments[4]=2 then curve at index 4 of dataset is assigned to cluster at index 2)

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
                        double dist = dF::discrete_frechet(this->centers[i], this->centers[j]);
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
        std::vector<curves::Curve2d> centers;          // cluster centers (constructor receives them initialized)
        std::vector<vector<curves::Curve2d>> clusters; // the clusters of curves

        CurveAssignor(int n_centers, std::vector<curves::Curve2d> centers, std::vector<curves::Curve2d> dataset) : n_centers(n_centers),
                                                                                                                   centers(centers),
                                                                                                                   dataset(dataset),
                                                                                                                   clusters(n_centers),
                                                                                                                   assignments_vec(dataset.size()) {}

        /* Assigns a nearest center to each point (part 1 of Lloyd's algorithm) using L2 distance.
         Stores both the current assignments and the current clusters in "assignments_vec" and "clusters" attributes respectively */
        void Lloyds()
        {
            int nearest_cntr;
            for (int i = 0; i < this->dataset.size(); ++i)
            {

                if (dataset[i].marked) // if item is marked, it is assigned in a cluster during reverse_assignement
                    continue;

                double min_d = dF::discrete_frechet(this->centers[0], this->dataset[i]);
                nearest_cntr = 0;

                for (int c = 1; c < centers.size(); ++c)
                {
                    double next_d = dF::discrete_frechet(this->centers[c], this->dataset[i]);
                    if (next_d < min_d)
                    {
                        min_d = next_d;
                        nearest_cntr = c;
                    }
                }

                this->assignments_vec[i] = nearest_cntr;
                this->clusters[nearest_cntr].push_back(this->dataset[i]); // push it into a cluster based on assigned center
                dataset[i].cluster = nearest_cntr;
            }
        }

        /* Assigns a nearest center to each point by reverse LSH range search and discrete Frechet distance.
         Expects a dFLSH::LSH object as parameter. */
        void Range_dfLSH(dFLSH::LSH &dflsh)
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

                            double dist_to_assigned = dF::discrete_frechet(this->centers[assigned_cluster], *item);

                            double dist_to_curr = dF::discrete_frechet(this->centers[c], *item);

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

            this->Lloyds();
        }
    };
}

#endif