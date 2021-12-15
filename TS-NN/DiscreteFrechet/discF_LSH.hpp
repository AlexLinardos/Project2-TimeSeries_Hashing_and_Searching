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

namespace dFLSH
{
    class LSH
    {
    private:
        std::default_random_engine eng;
        std::uniform_real_distribution<double> urd;
        vector<curves::Curve2d> dataset;
        int L;
        double delta;
        vector<curves::Curve2d *> **hashTables;
        int tableSize;

    public:
        LSH(vector<curves::Curve2d> &dataset, int L, double delta, int tableSize_divisor) : dataset(dataset),
                                                                                            L(L),
                                                                                            delta(delta),
                                                                                            tableSize(dataset.size() / tableSize_divisor),
                                                                                            eng(time(0) + clock()),
                                                                                            urd(0.0, delta){};

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

            // // apply padding if there were duplicates
            // double padding = 10000;
            // if (starting_size > clean_pis.size())
            // {
            //     cout << "I PADDED " << (starting_size - clean_pis.size()) << std::endl;
            //     for (int i = clean_pis.size(); i < starting_size; i++)
            //     {
            //         clean_pis.push_back(curves::Point2d(padding, padding));
            //     }
            // }

            return clean_pis;
        }

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

        void perform_hashing()
        {
        }
    };
}
#endif