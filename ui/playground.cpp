#include <iostream>
#include <string>
#include "Clustering_interface.hpp"

int main(int argc, char *argv[])
{
    Clusteri::Cluster_params params = Clusteri::Cluster_params(argc, argv);
    if (params.success == false)
    {
        return -1;
    }
    params.print_NN_params();

    return 0;
}