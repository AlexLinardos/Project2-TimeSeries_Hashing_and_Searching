#ifndef HASHING_HPP
#define HASHING_HPP
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <chrono>
#include "../utils.hpp"

using namespace std;

// hash function family H
class H
{
    int k; // number of H functions used in the algorithm
    int w; // window
    int d; // vector dimensions
    default_random_engine eng;
    normal_distribution<float> nd;
    uniform_real_distribution<double> urd;

public:
    H(int window, int dimensions, int k) : w(window), d(dimensions), k(k), eng(chrono::system_clock::now().time_since_epoch().count()), nd(0, 1), urd(0, w){};

    int produce_h(Item p)
    {
        vector<float> v;

        // create a vector v whose points follow the normal distribution
        for (int i = 0; i < d; i++)
        {
            v.push_back(nd(eng));
        }

        float t = urd(eng);

        return floor(((inner_product(p.xij.begin(), p.xij.end(), v.begin(), 0)) + t) / this->w);
    }

    vector<int> produce_k_h(Item p)
    {
        vector<int> h_vector;
        for (int i = 0; i < k; i++)
        {
            h_vector.push_back(produce_h(p));
        }
        return h_vector;
    }
};

#endif
