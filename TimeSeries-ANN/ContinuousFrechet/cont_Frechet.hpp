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
    // Curve vecToCurve(Item &item)
    // {
    //     Points points(2);
    //     for (int i = 0; i < item.xij.size(); i++)
    //     {
    //         // Coordinates coordinates;
    //         // coordinates.push_back((coordinate_t)item.xij[i]);
    //         // coordinates.push_back((coordinate_t)i);
    //         Point point(2);
    //         point.set(0, (coordinate_t)item.xij[i]);
    //         point.set(1, (coordinate_t)i);
    //         points.add(point);
    //     }

    //     // Curve * curve = new Curve(points, item.id);
    //     Curve curve(points, item.id);
    //     return curve;
    // }

    Curve fredCurve(curves::Curve2d &curve)
    {
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

        // Curve * fred_curve = new Curve(points, item.id);
        Curve fred_curve(points, curve.id);
        return fred_curve;
    }

    // double distance(Item &item1,  Item &item2)
    // {
    //     return Frechet::Continuous::distance(vecToCurve(item1), vecToCurve(item2)).value;
    // }

    double c_distance(curves::Curve2d &item1, curves::Curve2d &item2)
    {
        return Frechet::Continuous::distance(fredCurve(item1), fredCurve(item2)).value;
    }

    // searches for the exact nearest neighbour of the query curve using brute force approach
    std::pair<curves::Curve2d *, double> search_exactNN(curves::Curve2d &query, vector<curves::Curve2d> &dataset)
    {
        int size = query.data.size(); // length of curve

        // we will store current nearest neighbour in curr_NN along with its distance from query
        vector<double> dummy_vec;
        for (int i = 0; i < size; i++)
            dummy_vec.push_back(0.0);
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