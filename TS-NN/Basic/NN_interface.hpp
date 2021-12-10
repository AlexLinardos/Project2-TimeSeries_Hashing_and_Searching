#ifndef NN_INTERFACE_HPP
#define NN_INTERFACE_HPP
#include <iostream>
#include <string>

namespace NNi
{
    // class to hold parameters of Nearest Neighbor search
    // values are initialized to default
    class NN_params
    {
    private:
        std::string input_f;   // input file
        std::string query_f;   // query file
        int k;                 // number of hi hash functions for LSH or d' dimensions for HC
        int L;                 // number of hash tables for LSH
        int M;                 // max number of vectors to be checked for HC
        int probes;            // max number of probes to be checked for HC
        std::string output_f;  // ouput file
        std::string algorithm; // algorithm to be used for this run
        std::string metric;    // metric to be used for this run (only if algorithm = discrete or continuous)
        double delta;          // delta parameter for LSH of curves

    public:
        NN_params(int argc, char *argv[])
        {
            std::cout << "SUCCESS" << std::endl;
        }

        void print_NN_params()
        {
            std::cout << "[Current parameters for Neareset Neighbor search]" << std::endl
                      << "Input file: " << this->input_f << std::endl
                      << "Query file: " << this->query_f << std::endl
                      << "Output file: " << this->output_f << std::endl
                      << "Algorithm: " << this->algorithm << std::endl
                      << "Metric: " << this->metric << std::endl
                      << "k: " << this->k << std::endl
                      << "L: " << this->L << std::endl
                      << "M: " << this->M << std::endl
                      << "probes: " << this->probes << std::endl
                      << "delta: " << this->delta << std::endl;
        }
    };
}

#endif