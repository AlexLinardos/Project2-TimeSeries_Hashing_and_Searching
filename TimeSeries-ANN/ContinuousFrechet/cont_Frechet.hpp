#ifndef CON_FRECHET_HPP
#define CON_FRECHET_HPP
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../includes/utils.hpp"
#include "../../includes/curves.hpp"

#include "./Fred/frechet.hpp"

namespace cF
{
    Curve fredCurve(curves::Curve2d &curve) 
    { // converts a Curve2d object (our class) to a Curve object used in Fred folder from github to calculate continuous frechet
        Points points(1);
        for (int i = 0; i < curve.data.size(); i++)
        {
            Point point(1);
            point.set(0, (coordinate_t)curve.data[i].x);

            Point point2(1);
            point2.set(0, (coordinate_t)curve.data[i].y);

            points.add(point);
            points.add(point2);
        }

        Curve fred_curve(points, curve.id);
        return fred_curve;
    }

    double c_distance(curves::Curve2d &item1, curves::Curve2d &item2)
    { // "black box" function that converts our structures to Fred folder's structures to compute continuous frechet distance
        return Frechet::Continuous::distance(fredCurve(item1), fredCurve(item2)).value;
    }

    
    std::pair<curves::Curve2d *, double> search_exactNN(curves::Curve2d &query, vector<curves::Curve2d> &dataset)
    { // searches for the exact nearest neighbour of the query curve using brute force approach
        int size = query.data.size(); // length of curve

        // we will store current nearest neighbour in curr_NN along with its distance from query
        vector<double> dummy_vec;
        curves::Curve2d null_curve = curves::Curve2d("null", dummy_vec, dummy_vec);
        std::pair<curves::Curve2d *, double> curr_NN;
        curr_NN.first = &null_curve;
        curr_NN.second = std::numeric_limits<double>::max();

        // for each curve in the dataset
        for (int i = 0; i < dataset.size(); i++)
        {
            // calculate cont Frechet distance to it from given query
            double cfd = cF::c_distance(query, dataset[i]);
            // if nearer curve is found
            if (cfd < curr_NN.second)
            {
                // replace curr_NN
                curr_NN.first = &(dataset[i]);
                curr_NN.second = cfd;
            }
        }
        return curr_NN;
    }
}
#endif