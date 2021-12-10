#ifndef NN_INTERFACE_HPP
#define NN_INTERFACE_HPP
#include <iostream>
#include <string>
#include <map>

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
        // a map that will help us set the parameters from the command line
        std::map<std::string, std::string> param_set = {
            {"-i", "none"},
            {"-q", "none"},
            {"-k", "none"},
            {"-L", "none"},
            {"-M", "none"},
            {"-probes", "none"},
            {"-o", "none"},
            {"-algorithm", "none"},
            {"-metric", "none"},
            {"-delta", "none"}};
        std::map<std::string, std::string>::iterator it;

    public:
        NN_params(int argc, char *argv[])
        {
            for (int i = 1; i < argc; i = i + 2)
            {
                this->it = this->param_set.find(argv[i]); // try to find the parameter
                if (this->it == this->param_set.end())    // if user entered a parameter that does not exist
                {
                    std::cout << "Parameter " << argv[i] << " is not compatible. Please make sure you follow the format bellow: " << std::endl
                              << "./bin/search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes "
                              << "<int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete "
                              << "or continuous | only for -algorithm Frechet> -delta <double>" << std::endl;
                    break; // stop reading
                }
                else // if found
                {
                    this->it->second = argv[i + 1];
                    std::cout << "Parameter " << argv[i] << " set to " << argv[i + 1];
                }
            }
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