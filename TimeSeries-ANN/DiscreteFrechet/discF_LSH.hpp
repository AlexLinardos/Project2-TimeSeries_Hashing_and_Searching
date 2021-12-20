#ifndef DISCF_LSH_HPP
#define DISCF_LSH_HPP
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <math.h>
#include <set>
#include <limits>
#include "./disc_Frechet.hpp"
#include "../../includes/curves.hpp"
#include "../L2/LSH.hpp"

namespace dFLSH
{
    // creates associations between curve, grid-curve and vector x items
    class Association
    {
    public:
        curves::Curve2d *curve;
        vector<curves::Point2d> *grid_curve;
        vector<double> *x;

        Association(curves::Curve2d *curve, vector<curves::Point2d> *grid_curve, vector<double> *x) : curve(curve),
                                                                                                      grid_curve(grid_curve),
                                                                                                      x(x) {}
    };

    class LSH
    {
    private:
        std::default_random_engine eng;
        std::uniform_real_distribution<double> urd;
        vector<curves::Curve2d> *dataset;
        int L; // number of hash tables
        double delta;
        int tableSize; // size of each table
        int datasize;
        // int windowSize;
        vector<vector<curves::Point2d>> h_curves; // stores grid-curves
        vector<vector<double>> x_vecs;            // stores real vectors x
        std::vector<Association> **hashTables;    // Association* hashTables;
        vector<std::pair<double, double>> shifts; // stores shifts for all grids so we can snap queries on those grids
        double padding = 10000;                   // value that will replace duplicates
        G *g_family;                              // G hash family that is gonna be used for storing in 1d table

    public:
        LSH(vector<curves::Curve2d> *dataset, int L, double delta, int divisor_for_tableSize, bool querying_trick = false) : dataset(dataset),
                                                                                                                             L(L),
                                                                                                                             delta(delta),
                                                                                                                             shifts(L),
                                                                                                                             tableSize(dataset->size() / divisor_for_tableSize),
                                                                                                                             eng(time(0) + clock()),
                                                                                                                             urd(0.0, delta)
        {
            // // tune windowSize
            // std::random_device rd;                                          // only used once to initialise (seed) engine
            // std::mt19937 rng(rd());                                         // random-number engine used (Mersenne-Twister in this case)
            // std::uniform_int_distribution<int> uni(0, dataset->size() - 1); // guaranteed unbiased
            // int item_index_1;
            // int item_index_2;
            // double distance = 0;
            // /* For dataset->size()/4 samples we randomly choose two curves of the dataset and calculate their discrete frechet distance.
            //     We sum these disances and calculate the average. Then we multiply that average distance by a factor of our choice and we
            //     get the windowsize. That factor can be adjusted at search_main.cpp at the initialization of the LSH instance */
            // for (int i = 0; i < dataset->size() / 4; i++)
            // {
            //     item_index_1 = uni(rng);
            //     item_index_2 = uni(rng);
            //     while (item_index_1 == item_index_2)
            //         item_index_2 = uni(rng);
            //     distance += (dF::discrete_frechet((*dataset)[item_index_1], (*dataset)[item_index_2])) / (double)(dataset->size() / 4);
            // }
            // // this->windowSize = (int)floor(factor_for_windowSize * distance);
            // // std::cout << "windowSize: " << this->windowSize << std::endl;

            double windowSize = dF::mean_df_between_curves(*dataset);
            cout << "w " << windowSize << endl;

            this->g_family = new G(4, this->tableSize, windowSize, (*this->dataset)[0].data.size() * 2); // will be used for storing in 1d table

            // Initialize L hashTables, Grids(shifted_deltas)
            hashTables = new std::vector<Association> *[L];
            for (int i = 0; i < L; i++) // for every hashTable
            {
                hashTables[i] = new std::vector<Association>[tableSize];
                // shifted_deltas[i] =
            }
            // Hash all items in dataset and insert them into their buckets
            this->dataset_hashing();
            // for (int i = 0; i < this->L; i++)
            // {
            //     std::cout << "TABLE " << i << std::endl;
            //     for (int j = 0; j < this->tableSize; j++)
            //     {
            //         std::cout << "BUCKET " << j << " has " << this->hashTables[i][j].size() << std::endl;
            //     }
            // }
        };

        ~LSH()
        {
            delete this->g_family;
            for (int i = 0; i < this->L; i++)
            {
                delete[] hashTables[i];
            }
            delete[] hashTables;
        }

        // maps curve P to a grid
        vector<curves::Point2d> produce_h(curves::Curve2d curve, double tx, double ty)
        {
            int starting_size = curve.data.size();
            vector<curves::Point2d> snap_pis; // pi'

            // follow <<xi' = floor((x-t)/δ + 1/2)*δ + t>> formula to perform snapping
            for (int i = 0; i < starting_size; i++)
            {
                double snap_x = floor((curve.data[i].x - tx) / this->delta + 0.5) * this->delta + tx;
                double snap_y = floor((curve.data[i].y - ty) / this->delta + 0.5) * this->delta + ty;
                snap_pis.push_back(curves::Point2d(snap_x, snap_y));
            }

            // create new vector without consecutive duplicates
            vector<curves::Point2d> clean_pis;
            clean_pis.push_back(snap_pis[0]);
            for (int i = 1; i < starting_size; i++)
            {
                if ((snap_pis[i].x != clean_pis.back().x) || (snap_pis[i].y != clean_pis.back().y))
                {
                    clean_pis.push_back(snap_pis[i]);
                }
            }

            return clean_pis;
        }

        // maps curve's data to a grid
        vector<curves::Point2d> produce_h_for_data(vector<curves::Point2d> &curve, double tx, double ty)
        {
            int starting_size = curve.size();
            vector<curves::Point2d> snap_pis; // pi'

            // follow <<xi' = floor((x-t)/δ + 1/2)*δ + t>> formula to perform snapping
            for (int i = 0; i < starting_size; i++)
            {
                double snap_x = floor((curve[i].x - tx) / this->delta + 0.5) * this->delta + tx;
                double snap_y = floor((curve[i].y - ty) / this->delta + 0.5) * this->delta + ty;
                snap_pis.push_back(curves::Point2d(snap_x, snap_y));
            }

            // create new vector without consecutive duplicates
            vector<curves::Point2d> clean_pis;
            clean_pis.push_back(snap_pis[0]);
            for (int i = 1; i < starting_size; i++)
            {
                if ((snap_pis[i].x != clean_pis.back().x) || (snap_pis[i].y != clean_pis.back().y))
                {
                    clean_pis.push_back(snap_pis[i]);
                }
            }

            return clean_pis;
        }

        // concatenates points of grid-curve to produce a (real) vector x
        vector<double> concat_points(vector<curves::Point2d> points)
        {
            vector<double> result;
            for (int i = 0; i < points.size(); i++)
            {
                result.push_back(points[i].x);
                result.push_back(points[i].y);
            }
            return result;
        }

        // performs hashing to assing Association items to buckets
        void dataset_hashing()
        {
            // repeat L times (where is L is the number of tables)
            for (int i = 0; i < this->L; i++)
            {
                // produce 2 shift values (one for each dimension)
                double tx = urd(eng);
                double ty = urd(eng);
                // save them
                this->shifts[i] = std::make_pair(tx, ty);

                // for each curve
                for (int j = 0; j < this->dataset->size(); j++)
                {
                    int starting_size = (*dataset)[j].data.size();

                    // snap it to grid
                    this->h_curves.push_back(this->produce_h((*dataset)[j], tx, ty));
                    int new_size = this->h_curves.back().size();

                    // produce vector x
                    this->x_vecs.push_back(this->concat_points(this->h_curves.back()));

                    // apply padding if needed
                    if (starting_size > new_size)
                    {
                        for (int z = new_size; z < starting_size; z++)
                        {
                            // apply twice because curve is 2d
                            this->x_vecs.back().push_back(this->padding);
                            this->x_vecs.back().push_back(this->padding);
                        }
                    }

                    // create Association between curve, grid-curve and vector
                    Association ass = Association(&(*dataset)[j], &this->h_curves.back(), &this->x_vecs.back());
                    // create Item object so we can use produce_g from previous project
                    Item *item_for_g = new Item((*dataset)[j].id, this->x_vecs.back());
                    // get item hash value
                    unsigned int hval = (*this->g_family).produce_g(*item_for_g);
                    unsigned int pos = hval % (long unsigned)this->tableSize;
                    delete item_for_g; // we don't need it anymore
                    // store association in table
                    this->hashTables[i][pos].push_back(ass);
                }
            }
        }

        // searches for the approximate nearest neighbour of the query curve
        std::pair<curves::Curve2d *, double> search_ANN(curves::Curve2d &query, bool querying_trick, int threshold = 0)
        {
            int starting_size = query.data.size();

            // we will store current nearest neighbour in curr_NN along with its distance from query
            vector<double> dummy_vec;
            for (int i = 0; i < starting_size; i++)
                dummy_vec.push_back(0.0);
            curves::Curve2d null_curve = curves::Curve2d("null", dummy_vec, dummy_vec);
            std::pair<curves::Curve2d *, double> curr_NN;
            curr_NN.first = &null_curve;
            curr_NN.second = std::numeric_limits<double>::max();

            int searched = 0; // will be used to check if we reached threshold of checks
            // for each hash table
            for (int i = 0; i < this->L; i++)
            {
                // snap query it to grid
                vector<curves::Point2d> grid_curve = this->produce_h(query, this->shifts[i].first, this->shifts[i].second);
                int new_size = grid_curve.size();

                // produce vector x
                vector<double> x_vec = this->concat_points(grid_curve);

                // apply padding if needed
                if (starting_size > new_size)
                {
                    for (int z = new_size; z < starting_size; z++)
                    {
                        // apply twice because curve is 2d
                        x_vec.push_back(this->padding);
                        x_vec.push_back(this->padding);
                    }
                }

                // create Item object so we can use produce_g from previous project
                Item *item_for_g = new Item(query.id, x_vec);

                // find the bucket
                long unsigned id = this->g_family->produce_g(*item_for_g);
                long unsigned bucket = id % (long unsigned)this->tableSize;

                delete item_for_g; // we don't need it anymore

                // querying trick, check only curves with identical grid_curves
                int found = 0;
                if (querying_trick == true)
                {
                    // for each item in the bucket
                    for (int j = 0; j < this->hashTables[i][bucket].size(); j++)
                    {
                        if (identical_curves(h_curves.back(), *(this->hashTables[i][bucket][j].grid_curve)))
                        {
                            double dfd = dF::discrete_frechet(*(this->hashTables[i][bucket][j].curve), query);
                            // if nearer curve is found
                            if (dfd < curr_NN.second)
                            {
                                // replace curr_NN
                                curr_NN.first = this->hashTables[i][bucket][j].curve;
                                curr_NN.second = dfd;
                            }

                            cout << "QUERYING TRICK returned ann for " << query.id << endl;
                            found = 1;
                        }
                    }
                    if (found == 1)
                        return curr_NN;
                }

                // for each item in the bucket
                for (int j = 0; j < this->hashTables[i][bucket].size(); j++)
                {
                    // check if we bumped into same curve as current nearest before doing calculations
                    if (this->hashTables[i][bucket][j].curve->id != curr_NN.first->id)
                    {
                        double dfd = dF::discrete_frechet(*(this->hashTables[i][bucket][j].curve), query);
                        // if nearer curve is found
                        if (dfd < curr_NN.second)
                        {
                            // replace curr_NN
                            curr_NN.first = this->hashTables[i][bucket][j].curve;
                            curr_NN.second = dfd;
                        }
                        searched++;
                        if (threshold != 0 && searched >= threshold)
                            return curr_NN;
                    }
                }
            }
            // std::cout << "SEARCHED: " << searched << endl;
            return curr_NN;
        }

        // performs the LSH range search algorithm to find nearest neighbours of query in a given radius (using discrete Frechet distance)
        std::vector<std::pair<curves::Curve2d *, double>> RangeSearch(vector<curves::Point2d> &query, double radius, int threshold = 0)
        {
            std::vector<std::pair<curves::Curve2d *, double>> neighbours;

            int starting_size = query.size();

            int searched = 0; // will be used to check if we reached threshold of checks
            // for each hash table
            for (int i = 0; i < this->L; i++)
            {
                // snap query it to grid
                vector<curves::Point2d> grid_curve = this->produce_h_for_data(query, this->shifts[i].first, this->shifts[i].second);
                int new_size = grid_curve.size();

                // produce vector x
                vector<double> x_vec = this->concat_points(grid_curve);

                // apply padding if needed
                if (starting_size > new_size)
                {
                    for (int z = new_size; z < starting_size; z++)
                    {
                        // apply twice because curve is 2d
                        x_vec.push_back(this->padding);
                        x_vec.push_back(this->padding);
                    }
                }

                // create Item object so we can use produce_g from previous project
                Item *item_for_g = new Item("x_vec", x_vec);

                // find the bucket
                long unsigned id = this->g_family->produce_g(*item_for_g);
                long unsigned bucket = id % (long unsigned)this->tableSize;

                delete item_for_g; // we don't need it anymore

                // for each item in the bucket
                for (int j = 0; j < this->hashTables[i][bucket].size(); j++)
                {
                    // check if we bumped into same curve before doing calculations (and is already in neighbours vector)
                    for (int a = 0; a < neighbours.size(); a++)
                    {
                        /* In the "reverse assignment with range search using LSH" clustering algorithm we mark items when they are
                        assigned to a cluster so the next range search doesn't check them.*/
                        if ((this->hashTables[i][bucket][j].curve->id != neighbours[a].first->id) && (this->hashTables[i][bucket][j].curve->marked == false))
                        {
                            double dfd = dF::discrete_frechet_for_data(this->hashTables[i][bucket][j].curve->data, query);
                            // if curve is in radius
                            if (dfd < radius)
                            {
                                // add it to the neighbours of the query
                                neighbours.push_back(std::make_pair(hashTables[i][bucket][j].curve, dfd));
                            }
                            searched++;
                            if (threshold != 0 && searched >= threshold)
                                return neighbours;
                        }
                    }
                }
            }
            return neighbours;
        }
    };
}
#endif