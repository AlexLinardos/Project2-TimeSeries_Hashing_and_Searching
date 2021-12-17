#include <iostream>
#include <vector>
#include <string>
#include "../TimeSeries-ANN/DiscreteFrechet/disc_Frechet.hpp"
#include "../TimeSeries-ANN/ContinuousFrechet/contF_LSH.hpp"

// finds an optimal traversal between two curves using discrete Frechet distance
// optimal traversal will be returned in reverse so as to not waste time reversing it
vector<std::pair<int, int>> optimal_traversal(curves::Curve2d &p, curves::Curve2d &q)
{
    // get dynamic programming table from discrete Frechet
    double **c = dF::discrete_frechet(p, q);

    // initialise empty list(aka vector) of pairs
    std::vector<std::pair<int, int>> traversal;
    // pi = m1; qi = m2
    traversal.push_back(std::make_pair(p.data.size() - 1, q.data.size() - 1));
    int pi = traversal.back().first;
    int qi = traversal.back().second;
    while ((pi != 0) && (qi != 0))
    {
        // find next optimal step using the table of discrete Frechet
        int minIdx_p = pi - 1;
        int minIdx_q = qi;
        double mindist = c[pi - 1][qi];
        if (c[pi][qi - 1] < mindist)
        {
            minIdx_p = pi;
            minIdx_q = qi - 1;
            mindist = c[pi][qi - 1];
        }
        if (c[pi - 1][qi - 1] < mindist)
        {
            minIdx_p = pi - 1;
            minIdx_q = qi - 1;
            // mindist = c[pi - 1][qi - 1];
        }
        // add optimal step found to optimal traversal
        traversal.push_back(std::make_pair(minIdx_p, minIdx_q));
        pi = minIdx_p;
        qi = minIdx_q;
    }
    // corner-cases of endgame
    while (pi != 0) // if curve p still has steps to go but q has ended
    {
        traversal.push_back(std::make_pair(pi - 1, qi));
        pi--;
    }
    while (qi != 0) // if curve q still has steps to go but p has ended
    {
        traversal.push_back(std::make_pair(pi, qi - 1));
        qi--;
    }
    return traversal; // Warning: traversal is returned in reverse so as to save time
}

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
    cF::filter_curve(mean, 1.0);
    return mean;
}

// calculates mean curve of two given curves using discrete Frechet distance
std::vector<curves::Point2d> mean_of_curves(vector<curves::Curve2d> &curves)
{
    vector<vector<curves::Point2d>> curves_points;
    int i = 0;
    int remaining=curves.size();

    while(remaining>0)
    {
        if(remaining==1)
        {
            curves_points.push_back(curves[0].data);
            remaining -=1;
        }
        else
        {
            curves_points.push_back(mean_curve(curves[remaining-2], curves[remaining-1]));
            remaining -=2;
        }
    }
    //remaining=curves_points.size();

    while (curves_points.size()!=1)
    {
        cout << curves_points.size() << endl;
        remaining=curves_points.size();
        while(remaining>1)
        {
            cout << "remaining: " << remaining << endl;
            curves::Curve2d curve1 = curves::Curve2d("null", curves_points[remaining-2]);
            curves::Curve2d curve2 = curves::Curve2d("null", curves_points[remaining-1]);
            curves_points.push_back(mean_curve(curve1, curve2));
            curves_points.erase(curves_points.begin() + remaining-1);
            curves_points.erase(curves_points.begin() + remaining-2);
            remaining -=2;
        }
    }

    cout << curves_points.size() << endl;
    cout << curves_points[0].size() << endl;
    // for(int i=0; i<curves_points[0].size(); i++)
    // {
    //     cout << "[" << curves_points[0][i].x << ", " << curves_points[0][i].y << "] ";
    // }
    return curves_points[0];
}