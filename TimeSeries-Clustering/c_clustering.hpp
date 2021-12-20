#ifndef C_CLUSTERING_HPP
#define C_CLUSTERING_HPP

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
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
        Cli::Cluster_params params;
        int n_centers;  // number of centers to be initialized
        std::vector<curves::Curve2d> * dataset; // a dataset of curves

        vector<vector<curves::Point2d>> centers;
        vector<vector<curves::Curve2d>> clusters;
        std::vector<int> assignments_vec;     // shows the cluster to which each curve is assigned to
                                              // (example: if assignments[4]=2 then curve at index 4 of dataset is assigned to cluster at index 2)

        double mean_df;

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
        Clustering(Cli::Cluster_params &params, std::vector<curves::Curve2d> * dataset) : params(params),
                                                                                n_centers(params.clusters),
                                                                                clusters(params.clusters),
                                                                                assignments_vec(dataset->size()),
                                                                                dataset(dataset),
                                                                                mean_df(dF::mean_df_between_curves(*dataset)),
                                                                                eng(chrono::system_clock::now().time_since_epoch().count()),
                                                                                uid(0, dataset->size() - 1) {}

        void initialize_pp()
        {
            cout << "init " << endl;

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
            // cout << "assign " << clusters.size() << endl;
            int nearest_cntr;
            for (int i = 0; i < this->dataset->size(); ++i)
            {

                if ((*dataset)[i].marked) // if item is marked, it is assigned in a cluster during reverse_assignement
                    continue;

                double min_d = dF::discrete_frechet_for_data(this->centers[0], (*dataset)[i].data);
                nearest_cntr = 0;

                // cout << centers.size() << endl;

                for (int c = 1; c < centers.size(); ++c)
                {
                    // cout << c << endl;
                    double next_d = dF::discrete_frechet_for_data(this->centers[c], (*dataset)[i].data);
                    if (next_d < min_d)
                    {
                        min_d = next_d;
                        nearest_cntr = c;
                    }
                }

                // cout << "bef pushback " << nearest_cntr << endl;
                this->assignments_vec[i] = nearest_cntr;
                // cout << "bef pushback 2 " << assignments_vec[i] << endl;
                this->clusters[nearest_cntr].push_back((*dataset)[i]); // push it into a cluster based on assigned center
                // cout << "aft pushback " << endl;
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
            cout << "update " << endl;
            // we must calculate mean per cluster and make it the new center
            for (int i = 0; i < centers.size(); ++i)
            {
                centers[i] = mean_of_curves(clusters[i]);
            }
        }

        void perform_Lloyds(int max_iter)
        {
            cout << "perf lloyds " << endl;
            int iter = 1; // iterations
            this->initialize_pp();

            this->Lloyds_assignment();
            this->update_curve_centers();

            for (int i = 0; i < clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < clusters[i].size(); j++)
                {
                    std::cout << clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }

            vector<int> last_assignments(this->assignments_vec.size());
            do
            {
                last_assignments = this->assignments_vec;
                for (int i = 0; i < this->clusters.size(); ++i)
                {
                    this->clusters[i].clear(); // clear clusters so we can reassign the items
                }

                this->Lloyds_assignment();
                this->update_curve_centers();
                iter++;
                // iterate until assignments don't change or until we reach max_iter threshold
            } while ((!equal(assignments_vec.begin(), assignments_vec.end(), last_assignments.begin())) && iter < max_iter);

            std::cout << "AFTER LLOYDS " << std::endl;
            for (int i = 0; i < clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < clusters[i].size(); j++)
                {
                    std::cout << clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }

            cout << "............................................" << endl;
            cout << "Lloyd's algorithm ended after " << iter << " iterations" << endl;
        }

        void Reverse_Assignment_Clustering(int max_iterations)
        {
            double delta = delta_tuning(*dataset);

            // we must construct a discrete Frechet LSH object to pass as parameter in the assignment algorithm
            dFLSH::LSH *dflsh_object = new dFLSH::LSH(dataset, params.L, delta, 8);

            cout << "perf reverse " << endl;
            int iter = 1; // iterations

            this->initialize_pp();
            Range_dfLSH_assignment(*dflsh_object);
            update_curve_centers();
            vector<vector<curves::Point2d>> old_centers = this->centers;

            for (int i = 0; i < clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < clusters[i].size(); j++)
                {
                    std::cout << clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }

            do
            {
                for (int i = 0; i < this->clusters.size(); ++i)
                    this->clusters[i].clear(); // clear clusters so we can reassign the items
                for (int i = 0; i < dataset->size(); i++)
                    (*dataset)[i].marked = false;
                Range_dfLSH_assignment(*dflsh_object);
                old_centers = centers;
                update_curve_centers();
                iter++;

                cout << iter << " " << dF::max_centers_displacement(centers, old_centers) << endl;
            } while ((iter < max_iterations) && (dF::max_centers_displacement(centers, old_centers) > 1.0));
            
            for (int i = 0; i < clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < clusters[i].size(); j++)
                {
                    std::cout << clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }

            cout << "............................................" << endl;
            cout << "Reverse Assignment Cluestering ended after " << iter << " iterations" << endl;
        }

        // Silhouette of object at index i
        double silhouette(curves::Curve2d &curve)
        {
            vector<curves::Curve2d> * cluster = &(clusters[curve.cluster]); // find out at which cluster this item is assigned to

            //  calculate a(i) = average distance of i to objects in same cluster
            vector<double> distances;

            for (int j = 0; j < (*cluster).size(); ++j)
            {
                double dist = dF::discrete_frechet(curve, (*cluster)[j]);
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
                if (j != curve.cluster)
                {
                    double dist = dF::discrete_frechet_for_data(curve.data, this->centers[j]);
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
                double dist = dF::discrete_frechet(curve, (*cluster)[j]);
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