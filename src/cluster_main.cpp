#include <iostream>
#include <string>
#include "../ui/Clustering_interface.hpp"

int main(int argc, char *argv[])
{
    Cli::Cluster_params params = Cli::Cluster_params(argc, argv);
    if (params.success == false)
    {
    	std::cout<<"-----------------------------------------------------------"<<std::endl 
    		<< "Something went wrong while reading command line parameters." << std::endl
                  << "Please make sure you follow the format bellow: " << std::endl
                  << "./bin/cluster -i <input file> -c <configuration file> -o <output file> -update <Mean Frechet "
                  << "or Mean Vector> -assignment <Classic or LSH or Hypercube or LSH_Frechet> "
                  << "-complete <optional> -silhouette <optional>"<< std::endl;
        return -1;
    }
    params.print_NN_params();

    return 0;
}
