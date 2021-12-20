#ifndef VECTOR_CLUSTERING_HPP
#define VECTOR_CLUSTERING_HPP

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include "../includes/utils.hpp"
#include "../TimeSeries-ANN/L2/LSH.hpp"
#include "../TimeSeries-ANN/L2/HC.hpp"
#include "../ui/Clustering_interface.hpp"

namespace vector_cluster
{
    class Clustering
    {
        Cli::Cluster_params params;
        int n_centers;                         // number of centers to be initialized
        std::vector<Item> *dataset; // a dataset of curves
        int dimensions;  

        std::vector<int> assignments_vec;     // shows the cluster to which each curve is assigned to
                                              // (example: if assignments[4]=2 then curve at index 4 of dataset is assigned to cluster at index 2)

        std::default_random_engine eng; // an engine to help us pick first center randomnly
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
                double min_dist = EuclideanDistance(&centers[0], &(*dataset)[i], this->dimensions);
                if (min_dist > max_d)
                {
                    max_d = min_dist;
                }
                for (int cntr = 1; cntr < centers.size(); ++cntr)
                {
                    double dist = EuclideanDistance(&centers[cntr], &(*dataset)[i], this->dimensions);
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

        vector<Item> centers;
        vector<vector<Item>> clusters;

        Clustering(Cli::Cluster_params &params, std::vector<Item> * dataset) : params(params),
                                                                                dimensions(dataset[0].size()),
                                                                                n_centers(params.clusters),
                                                                                clusters(params.clusters),
                                                                                assignments_vec(dataset->size()),
                                                                                dataset(dataset),
                                                                                eng(chrono::system_clock::now().time_since_epoch().count()),
                                                                                uid(0, dataset->size() - 1) {}

        void initialize_pp()
        {
            // pick first centroid at random
            int rcentroid_index = this->uid(eng);
            centers.push_back((*dataset)[rcentroid_index]);

            // continue by picking more centroids until we reach the given number
            // centroids must be far away from each other
            while (centers.size() < this->n_centers)
            {
                double dist_sum = 0.0;
                std::vector<double> d(dataset->size(), 0.0); // distances of points to closest centroid
                std::vector<double> p(dataset->size(), 0.0); // probabilities of points to be chosen as next centroid

                this->calculate_min_dists(d, dist_sum); // calculate min distance of all points to the closest centroid for them

                this->calculate_probs(d, dist_sum, p); // calculate probability of each point to be chosen as next centroid
                double sum_prob = 0.0;
                for (int k = 0; k < p.size(); ++k)
                {
                    sum_prob += p[k];
                }

                std::discrete_distribution<> distribution{p.begin(), p.end()};
                int pick = distribution(eng);
                centers.push_back((*dataset)[pick]);
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

                double min_d = EuclideanDistance(&(centers[0]), &((*dataset)[i]), this->dimensions);
                nearest_cntr = 0;

                for (int c = 1; c < centers.size(); ++c)
                {
                    double next_d = EuclideanDistance(&centers[c], &((*dataset)[i]), this->dimensions);
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
        void Range_LSH_assignment(LSH &lsh)
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

            this->Lloyds_assignment();
        }

        void Range_HC_assignment(Hypercube &hc)
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

            this->Lloyds_assignment();
        }

        

        // update step of clustering for curves by calculation of mean curve
        void update_centers()
        {
            // we must calculate mean per cluster and make it the new center
            for (int i = 0; i < centers.size(); ++i)
            {
                vector<double> mean(this->dimensions, 0.0);
                int T = clusters[i].size(); // number of objects in cluster

                // calculate mean
                for (int j = 0; j < T; ++j)
                {
                    mean = vector_mean(mean, clusters[i][j].xij, dimensions, T);
                }

                centers[i].xij = mean;
            }
        }

        void Classic_Clustering(int max_iter)
        {
            std::cout << "[PERFORMING CLUSTERING] " << std::endl;
            int iter = 1; // iterations
            this->initialize_pp();

            this->Lloyds_assignment();
            this->update_centers();

            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            vector<int> last_assignments(this->assignments_vec.size());
            do
            {
                last_assignments = this->assignments_vec;
                for (int i = 0; i < this->clusters.size(); ++i)
                {
                    this->clusters[i].clear(); // clear clusters so we can reassign the items
                }

                this->Lloyds_assignment();
                this->update_centers();
                iter++;
                // iterate until assignments don't change or until we reach max_iter threshold
            } while ((!equal(assignments_vec.begin(), assignments_vec.end(), last_assignments.begin())) && iter < max_iter);

            // std::cout << "AFTER LLOYDS " << std::endl;
            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            cout << "............................................" << endl;
            cout << "Lloyd's algorithm ended after " << iter << " iterations" << endl;
        }

        void Reverse_Assignment_LSH_Clustering(int max_iterations)
        {
            // we must construct a discrete Frechet LSH object to pass as parameter in the assignment algorithm
            LSH_params lsh_params = LSH_params();
            lsh_params.k = this->params.k;
            lsh_params.L = this->params.L;
            LSH *lsh_object = new LSH(lsh_params, *dataset, 1.0, 8);

            int iter = 1; // iterations

            this->initialize_pp();
            Range_LSH_assignment(*lsh_object);
            update_centers();
            vector<Item> old_centers = this->centers;

            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            vector<int> last_assignments(this->assignments_vec.size());

            do
            {
                last_assignments = this->assignments_vec;
                for (int i = 0; i < this->clusters.size(); ++i)
                    this->clusters[i].clear(); // clear clusters so we can reassign the items
                for (int i = 0; i < dataset->size(); i++)
                    (*dataset)[i].marked = false;
                Range_LSH_assignment(*lsh_object);
                old_centers = centers;
                update_centers();
                iter++;
            } while ((!equal(assignments_vec.begin(), assignments_vec.end(), last_assignments.begin())) && iter < max_iterations);

            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            cout << "............................................" << endl;
            cout << "Reverse Assignment Cluestering ended after " << iter << " iterations" << endl;

            delete lsh_object;
        }

        void Reverse_Assignment_HC_Clustering(int max_iterations)
        {
            std::cout << "[PERFORMING CLUSTERING] " << std::endl;
            // we must construct a discrete Frechet LSH object to pass as parameter in the assignment algorithm
            Cube_params cube_params = Cube_params();
            cube_params.M = this->params.M;
            cube_params.k = this->params.k;
            cube_params.probes = this->params.probes;
            F f = F(cube_params.k);
            Hypercube * cube = new Hypercube(cube_params, (*dataset), 1.0, f.h_maps);

            int iter = 1; // iterations

            this->initialize_pp();
            Range_HC_assignment(*cube);
            update_centers();
            vector<Item> old_centers = this->centers;

            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            vector<int> last_assignments(this->assignments_vec.size());

            do
            {
                last_assignments = this->assignments_vec;
                for (int i = 0; i < this->clusters.size(); ++i)
                    this->clusters[i].clear(); // clear clusters so we can reassign the items
                for (int i = 0; i < dataset->size(); i++)
                    (*dataset)[i].marked = false;
                Range_HC_assignment(*cube);
                old_centers = centers;
                update_centers();
                iter++;
            } while ((!equal(assignments_vec.begin(), assignments_vec.end(), last_assignments.begin())) && iter < max_iterations);

            // for (int i = 0; i < clusters.size(); i++)
            // {
            //     std::cout << "CLUSTER " << i << std::endl;
            //     for (int j = 0; j < clusters[i].size(); j++)
            //     {
            //         std::cout << clusters[i][j].id << " ";
            //     }
            //     std::cout << std::endl;
            // }

            cout << "............................................" << endl;
            cout << "Reverse Assignment Cluestering ended after " << iter << " iterations" << endl;

            delete cube;
        }

        // Silhouette of object at index i
        double silhouette(Item &item)
        {
            vector<Item> *cluster = &clusters[item.cluster]; // find out at which cluster this item is assigned to

            //  calculate a(i) = average distance of i to objects in same cluster
            vector<double> distances;

            for (int j = 0; j < (*cluster).size(); ++j)
            {
                double dist = EuclideanDistance(&item, &(*cluster)[j], dimensions);
                distances.push_back(dist);
            }

            double a = 0.0;
            for (int j = 0; j < distances.size(); ++j)
            {
                a += distances[j] / (double)(distances.size());
            }

            // calculate b(i) = average distance of i to objects in next best (neighbor) cluster, i.e. cluster of 2nd closest centroid
            // first find the index of the next best cluster
            int best;
            double best_dist = std::numeric_limits<double>::max();
            for (int j = 0; j < this->centers.size(); ++j)
            {
                if (j != item.cluster)
                {
                    double dist = EuclideanDistance(&item, &this->centers[j], dimensions);
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
                double dist = EuclideanDistance(&item, &((*cluster)[j]), dimensions);
                distances.push_back(dist);
            }
            double b = 0.0;
            for (int j = 0; j < distances.size(); ++j)
            {
                b += distances[j] / (double)(distances.size());
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
                s = silhouette(clusters[c_index][i]);
                average += s / (double)n;
            }
            return average;
        }

        // metric to evaluate overall clustering
        double eval_clustering()
        {
            double average = 0.0;
            double s = 0.0;
            int n = this->dataset->size();
            for (int i = 0; i < n; ++i)
            {
                s = silhouette((*dataset)[i]);
                average += s / (double)n;
            }
            return s;
        }

    };
}
#endif