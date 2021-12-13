#include <iostream>
#include "NN_interface.hpp"
#include "LSH.hpp"

// ./bin/search -i <input file> -q <query file> -k 1 -L 2 -M 3 -probes 4 -o <output file> -algorithm <LSH> -metric <discrete> -delta 5.01
int main(int argc, char *argv[])
{
    NNi::NN_params params = NNi::NN_params(argc, argv);
    if (params.success == false)
    {
        std::cout << "Something went wrong while reading command line parameters." << std::endl
                  << "Please make sure you follow the format bellow: " << std::endl
                  << "./bin/search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes "
                  << "<int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete "
                  << "or continuous | only for -algorithm Frechet> -delta <double>" << std::endl;
        return -1;
    }
    params.print_NN_params();

    LSH_params lsh_params;
    lsh_params.input_file = params.input_f;
    lsh_params.query_file = params.query_f;
    lsh_params.k = params.k;
    lsh_params.L = params.L;
    lsh_params.out_file = params.output_f;
    /*std::cout << lsh_params.input_file << std::endl
              << lsh_params.query_file << std::endl
              << lsh_params.k << std::endl
              << lsh_params.L << std::endl
              << lsh_params.out_file << std::endl;*/
}
