#ifndef CURVES_HPP
#define CURVES_HPP
#include <vector>
#include <string>
#include <iostream>

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
        bool success = true;

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

#endif