#ifndef CONTF_LSH_HPP
#define CONTF_LSH_HPP
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <math.h>
#include <set>

#include "./cont_Frechet.hpp"
#include "../DiscreteFrechet/disc_Frechet.hpp"
#include "../../includes/curves.hpp"
#include "../L2/LSH.hpp"

namespace cF
{
    // minima maxima y removing from π˜ = v1, . . . , v` any vi s.t. vi ∈ [min{vi−1, vi+1}, max{vi−1, vi+1}]
    // στο min/max αφαιρεις ενα στοιχειο αν ειναι μεγαλυτερο απο το προηγουμενο του και μικροτερο απο το επομενο του ή αν ειναι ισο με ενα απο τα 2
    void minima_maxima(vector<double> &p)
    {
        int og_size = p.size();
        int i = 0;
        int cuts = 0;
        int total_cuts = 0;
        while (i + 2 < p.size())
        {
            cuts = 0;
            if (((p[i + 1] > min(p[i], p[i + 2])) && (p[i + 1] > max(p[i], p[i + 2]))) || p[i + 1] == p[i] || p[i + 1] == p[i + 2])
            {
                // cout << p[i] << " " << p[i+1] << " " << p[i+2] << endl;
                p.erase(p.begin() + i + 1); // erase element i+1
                cuts++;
                // cout << "cut1" << endl;
                total_cuts++;
            }
            if (cuts == 0)
                i++;
        }
        while (cuts != 0)
        {
            cuts = 0;
            i--;
            if (((p[i + 1] > min(p[i], p[i + 2])) && (p[i + 1] > max(p[i], p[i + 2]))) || p[i + 1] == p[i] || p[i + 1] == p[i + 2])
            {
                p.erase(p.begin() + i + 1); // erase element i+1
                cuts++;
                // cout << "cut2" << endl;
                total_cuts++;
            }
        }

        cout << "minima_maxima cut " << total_cuts << " points out of " << og_size << endl;
        return;
    }

    curves::Curve2d *filter_curve(curves::Curve2d &curve, double e) // for any consecutive points a, b, c, if |a − b| ≤ ε and |b − c| ≤ ε then remove b
    {
        curves::Curve2d *filtered_curve = new curves::Curve2d(curve.id);
        filtered_curve->data = curve.data;
        int i = 0;
        int cuts = 0;
        int total_cuts = 0;
        while (i + 2 < filtered_curve->data.size())
        {
            cuts = 0;
            if ((point2d_L2(filtered_curve->data[i], filtered_curve->data[i + 1]) <= e) && (point2d_L2(filtered_curve->data[i + 1], filtered_curve->data[i + 2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin() + i + 1); // erase element i+1
                cuts++;
                total_cuts++;
            }
            if (cuts == 0)
                i++;
        }
        while (cuts != 0)
        {
            cuts = 0;
            i--;
            if ((point2d_L2(filtered_curve->data[i], filtered_curve->data[i + 1]) <= e) && (point2d_L2(filtered_curve->data[i + 1], filtered_curve->data[i + 2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin() + i + 1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "filtering cut " << total_cuts << " points out of " << curve.data.size() << endl;
        return filtered_curve;
    }

    void filter_curve(vector<curves::Point2d> &curve, double e) // for any consecutive points a, b, c, if |a − b| ≤ ε and |b − c| ≤ ε then remove b
    {
        int og_size = curve.size();
        int i=0;
        int cuts=0;
        int total_cuts=0;
        while(i+2 < curve.size())
        {
            cuts=0;
            if((point2d_L2(curve[i], curve[i+1]) <= e) && (point2d_L2(curve[i+1], curve[i+2]) <= e))
            {
                curve.erase(curve.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
            if(cuts==0)
                i++;
        }
        while(cuts!=0)
        {
            cuts=0;
            i--;
            if((point2d_L2(curve[i], curve[i+1]) <= e) && (point2d_L2(curve[i+1], curve[i+2]) <= e))
            {
                curve.erase(curve.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "filtering cut " << total_cuts << " points out of "<< og_size << endl;
        return;
    }

    vector<curves::Curve2d> *filter_curves(vector<curves::Curve2d> &curve_dataset, double e)
    {
        vector<curves::Curve2d> *filtered_curve_dataset = new vector<curves::Curve2d>;
        for (int i = 0; i < curve_dataset.size(); i++)
        {
            cout << "curve#" << i << ": ";
            filtered_curve_dataset->push_back(*filter_curve(curve_dataset[i], e));
        }
        return filtered_curve_dataset;
    }
}

namespace cFLSH
{
    // creates associations between curve, grid-curve and vector x items
    class Association
    {
    public:
        curves::Curve2d *curve;
        vector<double> *grid_curve;
        vector<double> *x;

        Association(curves::Curve2d *curve, vector<double> *grid_curve, vector<double> *x) : curve(curve),
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
        int tableSize;
        int datasize;                      // size of each table
        vector<vector<double>> h_curves;       // stores grid-curves after snapping and minima_maxima
        vector<vector<double>> x_vecs;         // stores real vectors x (after padding)
        std::vector<Association> **hashTables; // Association* hashTables;
        double *t;                             // stores shifts for all grids
        double padding = 10000;
        G *g_family;

    public:
        LSH(vector<curves::Curve2d> *dataset, int L, double delta, int divisor_for_tableSize) : dataset(dataset),
                                                                                                                           L(L),
                                                                                                                           delta(delta),
                                                                                                                           tableSize(dataset->size() / divisor_for_tableSize),
                                                                                                                           eng(time(0) + clock()),
                                                                                                                           urd(0.0, delta)
        {
            // tune windowSize
            std::random_device rd;                                          // only used once to initialise (seed) engine
            std::mt19937 rng(rd());                                         // random-number engine used (Mersenne-Twister in this case)
            std::uniform_int_distribution<int> uni(0, dataset->size() - 1); // guaranteed unbiased
            int item_index_1;
            int item_index_2;
            double distance = 0;
            /* For dataset->size()/4 samples we randomly choose two points of the dataset and calculate their Euclidean distance.
                We sum these disances and calculate the average. Then we multiply that average distance by a factor of our choice and we
                get the windowsize. That factor can be adjusted at LSH.cpp at the initialization of the LSH instance */
            for (int i = 0; i < dataset->size() / 4; i++)
            {
                item_index_1 = uni(rng);
                item_index_2 = uni(rng);
                while (item_index_1 == item_index_2)
                    item_index_2 = uni(rng);
                distance += dF::discrete_frechet((*dataset)[item_index_1], (*dataset)[item_index_2]) / (double)(dataset->size() / 4);
                // distance += (cF::c_distance((*dataset)[item_index_1], (*dataset)[item_index_2])) / (double)(dataset->size() / 4);
            }

            this->g_family = new G(4, this->tableSize, distance, (*this->dataset)[0].data.size() * 2); // will be used for storing in 1d table
            
            // Initialize L hashTables, Grids(shifted_deltas)
            hashTables = new std::vector<Association> *[L];
            t = new double[L];
            for (int i = 0; i < L; i++) // for every hashTable
            {
                hashTables[i] = new std::vector<Association>[tableSize];
                t[i] = urd(eng);
            }

            // Hash all items in dataset and insert them into their buckets

            this->dataset_hashing();
            for (int i = 0; i < this->L; i++)
            {
                std::cout << "TABLE " << i << std::endl;
                for (int j = 0; j < this->tableSize; j++)
                {
                    std::cout << "BUCKET " << j << " has " << this->hashTables[i][j].size() << std::endl;
                }
            }
        };

        ~LSH()
        {
            delete this->g_family;
            for (int i = 0; i < this->L; i++)
            {
                delete[] hashTables[i];
            }
            delete[] hashTables;
            delete t;
        }

        // maps curve P to a 1D grid
        vector<double> snap_to_1dgrid(curves::Curve2d &curve, double t)
        {
            int starting_size = curve.data.size();
            vector<double> snap_pis;

            // follow <<xi' = floor((x-t)/δ + 1/2)*δ + t>> formula to perform snapping
            for (int i = 0; i < starting_size; i++)
            {
                double snap_x = floor((curve.data[i].x - t) / this->delta + 0.5) * this->delta + t;
                double snap_y = floor((curve.data[i].y - t) / this->delta + 0.5) * this->delta + t;
                snap_pis.push_back(snap_x);
                snap_pis.push_back(snap_y);
            }

            // create new vector without consecutive duplicates
            vector<double> clean_pis;
            clean_pis.push_back(snap_pis[0]);
            for (int i = 1; i < starting_size; i++)
            {
                if ((snap_pis[i] != clean_pis.back()))
                {
                    clean_pis.push_back(snap_pis[i]);
                }
            }

            return clean_pis;
        }

        // minima maxima y removing from π˜ = v1, . . . , v` any vi s.t. vi ∈ [min{vi−1, vi+1}, max{vi−1, vi+1}]
        // minima_maxima(vector<double> &p)

        // performs hashing to assing Association items to buckets
        void dataset_hashing()
        {
            // repeat L times (where is L is the number of tables)
            for (int i = 0; i < this->L; i++)
            {
                // for each curve
                for (int j = 0; j < this->dataset->size(); j++)
                {
                    int starting_size = (*dataset)[j].data.size();

                    // snap it to grid (remove consecutive duplicates and minima_maxima)
                    this->h_curves.push_back(this->snap_to_1dgrid((*dataset)[j], t[i]));
                    cF::minima_maxima(h_curves.back());

                    // produce vector x
                    this->x_vecs.push_back(h_curves.back());

                    // apply padding if needed
                    int new_size = this->x_vecs.back().size();
                    if (starting_size > new_size)
                    {
                        for (int z = new_size; z < starting_size; z++)
                        {
                            this->x_vecs.back().push_back(this->padding);
                        }
                    }

                    // create Association between curve, grid-curve and x_vector
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
                // snap it to grid (remove consecutive duplicates and minima_maxima)
                this->h_curves.push_back(this->snap_to_1dgrid(query, t[i]));
                cF::minima_maxima(h_curves.back());

                // produce vector x
                this->x_vecs.push_back(h_curves.back());

                // apply padding if needed
                int new_size = this->x_vecs.back().size();
                if (starting_size > new_size)
                {
                    for (int z = new_size; z < starting_size; z++)
                    {
                        this->x_vecs.back().push_back(this->padding);
                    }
                }

                // create Item object so we can use produce_g from previous project
                Item *item_for_g = new Item(query.id, this->x_vecs.back());

                // find the bucket
                long unsigned hval = (*this->g_family).produce_g(*item_for_g);
                long unsigned bucket = hval % (long unsigned)this->tableSize;
                delete item_for_g; // we don't need it anymore

                // querying trick, if curve with identical grid_curve exists in bucket, return it as the ann
                if(querying_trick==true)
                {
                    // for each item in the bucket
                    for (int j = 0; j < this->hashTables[i][bucket].size(); j++)
                    {
                        if (identical_curves(h_curves.back(), *(this->hashTables[i][bucket][j].grid_curve)))
                        {
                            // replace curr_NN
                            curr_NN.first = this->hashTables[i][bucket][j].curve;
                            curr_NN.second = cF::c_distance(*(this->hashTables[i][bucket][j].curve), query);
                            cout << "QUERYING TRICK returned ann for " << query.id << endl;
                            return curr_NN;
                        }
                    }
                }

                // for each item in the bucket
                for (int j = 0; j < this->hashTables[i][bucket].size(); j++)
                {
                    // check if we bumped into same curve as current nearest before doing calculations
                    if (this->hashTables[i][bucket][j].curve->id != curr_NN.first->id)
                    {
                        double cfd = cF::c_distance(*(this->hashTables[i][bucket][j].curve), query);
                        // if nearer curve is found
                        if (cfd < curr_NN.second)
                        {
                            // replace curr_NN
                            curr_NN.first = this->hashTables[i][bucket][j].curve;
                            curr_NN.second = cfd;
                        }
                        searched++;
                        if (threshold != 0 && searched >= threshold)
                            return curr_NN;
                    }
                }
            }
            std::cout << "SEARCHED: " << searched << endl;
            return curr_NN;
        }
    };
}

#endif