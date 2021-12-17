#ifndef CURVES_HPP
#define CURVES_HPP
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "./optimal_traversal.hpp"

namespace curves
{
    // creates a point for a 2-dimensional space
    class Point2d
    {
    public:
        double x;
        double y;

        Point2d(double x, double y) : x(x), y(y){};
    };

    // calculates euclidean distance between two Point2d objects
    double point2d_L2(Point2d &point1, Point2d &point2)
    {
        double dist = pow((point1.x - point2.x), 2) + pow((point1.y - point2.y), 2);
        return sqrt(dist);
    }

    // creates a curve for a 2-dimensional space
    class Curve2d
    {
    public:
        std::string id;
        std::vector<Point2d> data;
        bool success = true;

        Curve2d(std::string id) : id(id) {}

        Curve2d(std::string id, std::vector<double> vector1, std::vector<double> vector2) : id(id)
        {

            if (vector1.size() != vector2.size())
            {
                std::cout << "Please provider data of same dimensions in order to construct curve." << std::endl;
                this->success = false;
            }
            else
            {
                for (int i = 0; i < vector1.size(); i++)
                {
                    this->data.push_back(Point2d(vector1[i], vector2[i]));
                }
            }
        }
    };

    // calculates mean curve of two given curves using discrete Frechet distance
    std::vector<curves::Point2d> mean_curve(curves::Curve2d &p, curves::Curve2d &q)
    {
        std::vector<curves::Point2d> mean;
        std::vector<std::pair<int, int>> opt_traversal = optimal_traversal(p, q);

        // iterate in reverse because optimal traversal will be return in reverse from optimal_traversal(p,q)
        for (int t = opt_traversal.size() - 1; t >= 0; t--)
        {
            std::pair<int, int> traversal_indexes = opt_traversal[t];
            double pit_x = p.data[traversal_indexes.first].x;
            double pit_y = p.data[traversal_indexes.first].y;
            double qit_x = q.data[traversal_indexes.second].x;
            double qit_y = q.data[traversal_indexes.second].y;
            mean.push_back(curves::Point2d((pit_x + qit_x) / 2, (pit_y + qit_y) / 2));
        }
        return mean;
    }
}

#endif