#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include "../includes/curves.hpp"
#include "../includes/utils.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/disc_Frechet.hpp"

namespace init
{
    // helps us implement initialization++ for a dataset of vectors
    class VectorInitializer
    {
    private:
        int n_centers;                  // number of centers to be initialized
        std::vector<Item> dataset;      // a dataset of vectors
        int dimensions;                 // dimensions of vectors in the dataset
        std::default_random_engine eng; // an engine to help us pick first center randomnly
        std::uniform_int_distribution<int> uid;
        std::vector<Item> centers; // cluster centers initialized with initialization++

        /* Calculates the L2 distance of each point in dataset from its nearest centroid
            Parameters are an std::vector containing the current centroids, and std::vector in which the min distance from nearest centroid
            of each corresponding dataset item will be stored and a sum accumulator "dist_sum" of all those distances (needed for calculate_probs).*/
        void calculate_min_dists(vector<Item> &centroids, vector<double> &d, double &dist_sum)
        {
            double max_d = 0.0; // for normalization
            // for each point in the dataset
            for (int i = 0; i < this->dataset.size(); ++i)
            {
                // calculate L2 distance from nearest centoid
                double min_dist = EuclideanDistance(&centroids[0], &this->dataset[i], this->dimensions);
                if (min_dist > max_d)
                {
                    max_d = min_dist;
                }
                for (int cntr = 1; cntr < centroids.size(); ++cntr)
                {
                    double dist = EuclideanDistance(&centroids[cntr], &this->dataset[i], this->dimensions);
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

    public:
        VectorInitializer(int n_centers, std::vector<Item> dataset, int dimensions) : n_centers(n_centers),
                                                                                      dataset(dataset),
                                                                                      dimensions(dimensions),
                                                                                      eng(chrono::system_clock::now().time_since_epoch().count()),
                                                                                      uid(0, dataset.size() - 1) {}

        // initialization ++ of cluster centroids (final centroids are both stored as attribute and returned)
        std::vector<Item> initialize_pp()
        {
            std::vector<Item> centroids;

            // pick first centroid at random
            int rcentroid_index = this->uid(eng);
            centroids.push_back(dataset[rcentroid_index]);

            // continue by picking more centroids until we reach the given number
            // centroids must be far away from each other
            while (centroids.size() < this->n_centers)
            {
                double dist_sum = 0.0;
                std::vector<double> d(dataset.size(), 0.0); // distances of points to closest centroid
                std::vector<double> p(dataset.size(), 0.0); // probabilities of points to be chosen as next centroid

                this->calculate_min_dists(centroids, d, dist_sum); // calculate min distance of all points to the closest centroid for them
                // cout << "SUM (after normalization) = " << dist_sum << endl;

                this->calculate_probs(d, dist_sum, p); // calculate probability of each point to be chosen as next centroid
                double sum_prob = 0.0;
                for (int k = 0; k < p.size(); ++k)
                {
                    sum_prob += p[k];
                }

                std::discrete_distribution<> distribution{p.begin(), p.end()};
                int pick = distribution(eng);
                centroids.push_back(dataset[pick]);
            }
            this->centers = centroids;
            return centroids;
        }
    };

    // helps us implement initialization++ for a dataset of curves
    class CurveInitializer
    {
    private:
        int n_centers;                        // number of centers to be initialized
        std::vector<curves::Curve2d> dataset; // a dataset of curves
        std::default_random_engine eng;       // an engine to help us pick first center randomnly
        std::uniform_int_distribution<int> uid;
        std::vector<curves::Curve2d> centers; // cluster centers initialized with initialization++

        /* Calculates the Frechet distance of each point in dataset from its nearest centroid
            Parameters are an std::vector containing the current centroids, and std::vector in which the min distance from nearest centroid
            of each corresponding dataset item will be stored and a sum accumulator "dist_sum" of all those distances (needed for calculate_probs).*/
        void calculate_min_dists(vector<curves::Curve2d> &centroids, vector<double> &d, double &dist_sum)
        {
            double max_d = 0.0; // for normalization
            // for each point in the dataset
            for (int i = 0; i < this->dataset.size(); ++i)
            {
                // calculate Frechet distance from nearest centoid
                double min_dist = dF::discrete_frechet(centroids[0], this->dataset[i]);
                if (min_dist > max_d)
                {
                    max_d = min_dist;
                }
                for (int cntr = 1; cntr < centroids.size(); ++cntr)
                {
                    double dist = dF::discrete_frechet(centroids[cntr], this->dataset[i]);
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

    public:
        CurveInitializer(int n_centers, std::vector<curves::Curve2d> dataset) : n_centers(n_centers),
                                                                                dataset(dataset),
                                                                                eng(chrono::system_clock::now().time_since_epoch().count()),
                                                                                uid(0, dataset.size() - 1) {}

        // initialization ++ of cluster centroids (final centroids are both stored as attribute and returned)
        std::vector<curves::Curve2d> initialize_pp()
        {
            std::vector<curves::Curve2d> centroids;

            // pick first centroid at random
            int rcentroid_index = this->uid(eng);
            centroids.push_back(dataset[rcentroid_index]);

            // continue by picking more centroids until we reach the given number
            // centroids must be far away from each other
            while (centroids.size() < this->n_centers)
            {
                double dist_sum = 0.0;
                std::vector<double> d(dataset.size(), 0.0); // distances of points to closest centroid
                std::vector<double> p(dataset.size(), 0.0); // probabilities of points to be chosen as next centroid

                this->calculate_min_dists(centroids, d, dist_sum); // calculate min distance of all points to the closest centroid for them
                // cout << "SUM (after normalization) = " << dist_sum << endl;

                this->calculate_probs(d, dist_sum, p); // calculate probability of each point to be chosen as next centroid
                double sum_prob = 0.0;
                for (int k = 0; k < p.size(); ++k)
                {
                    sum_prob += p[k];
                }

                std::discrete_distribution<> distribution{p.begin(), p.end()};
                int pick = distribution(eng);
                centroids.push_back(dataset[pick]);
            }
            this->centers = centroids;
            return centroids;
        }
    };
}

#endif