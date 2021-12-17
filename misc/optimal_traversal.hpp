#ifndef OPT_TRAVERSAL_HPP
#define OPT_TRAVERSAL_HPP
#include <iostream>
#include <vector>
#include <string>
#include "./curves.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/disc_Frechet.hpp"

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

#endif