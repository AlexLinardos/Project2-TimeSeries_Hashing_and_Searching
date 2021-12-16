#ifndef DISCF_LSH_HPP
#define DISCF_LSH_HPP
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <math.h>
#include <set>
#include "disc_Frechet.hpp"
#include "../curves.hpp"
#include "../Basic/LSH.hpp"

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
        vector<curves::Curve2d> dataset;
        int L; // number of hash tables
        double delta;
        int tableSize;                            // size of each table
        vector<vector<curves::Point2d>> h_curves; // stores grid-curves
        vector<vector<double>> x_vecs;            // stores real vectors x

    public:
        Association **hashTables;
        LSH(vector<curves::Curve2d> &dataset, int L, double delta, int tableSize_divisor) : dataset(dataset),
                                                                                            L(L),
                                                                                            delta(delta),
                                                                                            tableSize(dataset.size() / tableSize_divisor),
                                                                                            eng(time(0) + clock()),
                                                                                            urd(0.0, delta)
        {
            this->hashTables = new Association *[L];
            for (int i = 0; i < L; i++)
            {
                this->hashTables[i] = new Association[tableSize];
            }
        };
        // maps curve P to a grid
        vector<curves::Point2d> produce_h(curves::Curve2d curve)
        {
            int starting_size = curve.data.size();
            vector<curves::Point2d> snap_pis; // pi'
            // produce 2 shift values (one for each dimension)
            double tx = urd(eng);
            double ty = urd(eng);

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

        ~LSH()
        {
            for (int i = 0; i < this->L; i++)
            {
                delete[] hashTables[i];
            }
            delete[] hashTables;
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
            double padding = 10000;

            // repeat L times (where is L is the number of tables)
            for (int i = 0; i < this->L; i++)
            {
                // will be used for storing in 1d table
                G g_family = G(4, this->tableSize, 600); // w=600

                // for each curve
                for (int j = 0; j < this->dataset.size(); j++)
                {
                    int starting_size = this->dataset[j].data.size();
                    // snap it to grid
                    this->h_curves.push_back(this->produce_h(this->dataset[j]));
                    int new_size = this->h_curves.back().size();
                    // apply padding if needed
                    if (starting_size > new_size)
                    {
                        for (int z = new_size; z < starting_size; z++)
                        {
                            this->h_curves.back().push_back(curves::Point2d(padding, padding));
                        }
                    }
                    // produce vector x
                    this->x_vecs.push_back(this->concat_points(this->h_curves.back()));
                    // create Association between curve, grid-curve and vector
                    Association ass = Association(&this->dataset[j], &this->h_curves.back(), &this->x_vecs.back());
                    // create Item object so we can use produce_g from previous project
                    Item item_for_g = Item(this->dataset[j].id, this->x_vecs.back());
                    // get item hash value
                    unsigned int hval = g_family.produce_g(item_for_g);
                    // store it in table
                    this->hashTables[L][hval] = ass;
                }
            }
        }
    };
}
#endif