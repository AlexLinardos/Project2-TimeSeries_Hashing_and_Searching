#ifndef DISCF_LSH_HPP
#define DISCF_LSH_HPP
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include "disc_Frechet.hpp"

namespace dFLSH
{
    class Grid
    {
    private:
        std::default_random_engine eng;
        std::uniform_real_distribution<double> urd;

    public:
        double delta;
        Grid(double delta) : delta(delta), eng(time(0) + clock()), urd(0.0, delta){};
    };
}

#endif