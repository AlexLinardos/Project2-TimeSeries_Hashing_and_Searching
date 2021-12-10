#include <iostream>
#include "NN_interface.hpp"

int main(int argc, char *argv[])
{
    NNi::NN_params params = NNi::NN_params(argc, argv);
    params.print_NN_params();
}