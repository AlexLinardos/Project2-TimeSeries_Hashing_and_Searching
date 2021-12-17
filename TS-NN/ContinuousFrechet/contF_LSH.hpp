#ifndef CONTF_LSH_HPP
#define CONTF_LSH_HPP
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <math.h>
#include <set>

#include "cont_Frechet.hpp"
#include "../DiscreteFrechet/disc_Frechet.hpp"
#include "../curves.hpp"
#include "../Basic/LSH.hpp"

namespace cF
{
    // minima maxima y removing from π˜ = v1, . . . , v` any vi s.t. vi ∈ [min{vi−1, vi+1}, max{vi−1, vi+1}] 
    // στο min/max αφαιρεις ενα στοιχειο αν ειναι μεγαλυτερο απο το προηγουμενο του και μικροτερο απο το επομενο του ή αν ειναι ισο με ενα απο τα 2
    void minima_maxima(vector<double> &p)
    {
        int og_size = p.size();
        int i=0;
        int cuts=0;
        int total_cuts=0;
        while(i+2 < p.size())
        {
            cuts=0;
            if ((p[i+1] >= min(p[i], p[i+2])) || (p[i+1] >= max(p[i], p[i+2])))
            {
                p.erase(p.begin()+i+1); // erase element i+1
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
            if ((p[i+1] >= min(p[i], p[i+2])) || (p[i+1] >= max(p[i], p[i+2])))
            {
                p.erase(p.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "minima_maxima cut " << total_cuts << " points out of "<< og_size << endl;
        return;
    }

    
    curves::Curve2d * filter_curve(curves::Curve2d &curve, double e) // for any consecutive points a, b, c, if |a − b| ≤ ε and |b − c| ≤ ε then remove b
    {
        curves::Curve2d * filtered_curve = new curves::Curve2d(curve.id);
        filtered_curve->data=curve.data;
        int i=0;
        int cuts=0;
        int total_cuts=0;
        while(i+2 < filtered_curve->data.size())
        {
            cuts=0;
            if((point2d_L2(filtered_curve->data[i], filtered_curve->data[i+1]) <= e) && (point2d_L2(filtered_curve->data[i+1], filtered_curve->data[i+2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin()+i+1); // erase element i+1
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
            if((point2d_L2(filtered_curve->data[i], filtered_curve->data[i+1]) <= e) && (point2d_L2(filtered_curve->data[i+1], filtered_curve->data[i+2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "filtering cut " << total_cuts << " points out of "<< curve.data.size() << endl;
        return filtered_curve;
    }

    vector<curves::Curve2d> * filter_curves(vector<curves::Curve2d> &curve_dataset, double e)
    {
        vector<curves::Curve2d> * filtered_curve_dataset = new vector<curves::Curve2d>;
        for(int i=0; i < curve_dataset.size(); i++)
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
        int tableSize;
        int datasize;
        int windowSize;                           // size of each table
        vector<vector<curves::Point2d>> h_curves; // stores grid-curves
        vector<vector<double>> x_vecs;            // stores real vectors x
        std::vector<Association> **hashTables;    // Association* hashTables;
        double *t;                  // stores shifts for all grids

    public:
        LSH(vector<curves::Curve2d> *dataset, int L, double delta, int factor_for_windowSize, int divisor_for_tableSize) : dataset(dataset),
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
                distance += (dF::discrete_frechet((*dataset)[item_index_1], (*dataset)[item_index_2])[dataset->size() - 1][dataset->size() - 1]) / (double)(dataset->size() / 4);
                //distance += (cF::c_distance((*dataset)[item_index_1], (*dataset)[item_index_2])) / (double)(dataset->size() / 4);
            }
            windowSize = factor_for_windowSize * (int)distance;
            std::cout << "windowSize :" << windowSize << std::endl;

            // Initialize L hashTables, Grids(shifted_deltas)
            hashTables = new std::vector<Association> *[L];
            t = new double[L];
            for (int i = 0; i < L; i++) // for every hashTable
            {
                hashTables[i] = new std::vector<Association>[tableSize];
                // t[i] =
            }

            // Hash all items in dataset and insert them into their buckets


            // this->dataset_hashing();
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
            for (int i = 0; i < this->L; i++)
            {
                delete[] hashTables[i];
            }
            delete[] hashTables;
            delete t;
        }

        // maps curve P to a 1D grid
        vector<double> snap_to_1dgrid(curves::Curve2d * curve, double t)
        {
            int starting_size = curve->data.size();
            vector<double> snap_pis;

            // follow <<xi' = floor((x-t)/δ + 1/2)*δ + t>> formula to perform snapping
            for (int i = 0; i < starting_size; i++)
            {
                double snap_x = floor((curve->data[i].x - t) / this->delta + 0.5) * this->delta + t;
                double snap_y = floor((curve->data[i].y - t) / this->delta + 0.5) * this->delta + t;
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
    };
}






# endif