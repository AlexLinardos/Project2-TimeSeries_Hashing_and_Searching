#ifndef CURVES_HPP
#define CURVES_HPP
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <random>
#include <time.h>
#include "./utils.hpp"

namespace curves
{
    class Point2d
    {
    public:
        double x;
        double y;

        Point2d(double x, double y) : x(x), y(y){};
    };

    class Curve2d
    {
    public:
        std::string id;
        std::vector<Point2d> data;
        bool success = true; // checks that object was constructed succesfully using 2 vectors of same length

        int cluster = 0; // index of cluster to which this item is assigned
        // int cluster2 = -1; // second nearest cluster

        bool null = false;    // flag if is NULL item created to initialize a pair vector for knn and brute foce algorithms
        bool claimed = false; // will be used in reverse assignment to indicate if item has been claimed by a ball (in order to resolve conflicts)
        bool marked = false;  // will be used in reverse assignment to indicate item has been assigned to a cluster

        Curve2d(std::string id) : id(id) {}

        Curve2d(std::string id, std::vector<Point2d> &p_data) : id(id)
        {
            for (int i = 0; i < p_data.size(); i++)
            {
                this->data.push_back(p_data[i]);
            }
        }

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
}

double point2d_L2(curves::Point2d &point1, curves::Point2d &point2)
{
    double dist = pow((point1.x - point2.x), 2) + pow((point1.y - point2.y), 2);
    return sqrt(dist);
}

double avg_p_dist(curves::Curve2d &curve) // calculates the average distance between the vertices of a curve
{
    double avg = 0.0;
    for (int i = 0; i < curve.data.size() - 1; i++)
    {
        avg += point2d_L2(curve.data[i], curve.data[i + 1]) / curve.data.size();
    }
    return avg;
}

double delta_tuning(std::vector<curves::Curve2d> &curves) // calculates the average distance between the vertices of a curve
{
    std::random_device rd;                                        // only used once to initialise (seed) engine
    std::mt19937 rng(rd());                                       // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0, curves.size() - 1); // guaranteed unbiased

    int curve_index_1;
    int curve_index_2;

    double avg = 0.0;

    /*
        For dataset.size()/4 samples we randomly choose a curve and calculate the average distance between its' vertices.
        We sum these average disances and calculate the average.
    */

    for (int i = 0; i < curves.size() / 4; i++)
    {
        curve_index_1 = uni(rng);
        // curve_index_1 = uni(rng);
        // while (curve_index_1 == curve_index_1)
        //     curve_index_1 = uni(rng);
        avg += avg_p_dist(curves[curve_index_1]) / (curves.size() / 4);
    }
    std::cout << "delta: " << avg << std::endl;
    return avg;
}

bool identical_curves(vector<curves::Point2d> &curve1, vector<curves::Point2d> &curve2)
{
    // int curve_size = curve1.data.size();

    if (curve1.size() != curve2.size())
        return false;

    for (int i = 0; i < curve1.size(); i++)
    {
        if ((curve1[i].x != curve2[i].x) || (curve1[i].y != curve2[i].y))
            return false;
    }
    return true;
}

bool identical_curves(vector<double> &curve1, vector<double> &curve2)
{
    // int curve_size = curve1.size();

    if (curve1.size() != curve2.size())
        return false;

    for (int i = 0; i < curve1.size(); i++)
    {
        if ((curve1[i] != curve2[i]))
            return false;
    }
    return true;
}

#endif