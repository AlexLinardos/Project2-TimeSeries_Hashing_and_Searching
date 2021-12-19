#include <iostream>
#include <string>
#include <vector>
#include "../ui/Clustering_interface.hpp"
#include "../TimeSeries-Clustering/initialization.hpp"

int main(int argc, char *argv[])
{
    Cli::Cluster_params params = Cli::Cluster_params(argc, argv);
    if (params.success == false)
    {
        std::cout << "-----------------------------------------------------------" << std::endl
                  << "Something went wrong while reading command line parameters." << std::endl
                  << "Please make sure you follow the format bellow: " << std::endl
                  << "./bin/cluster -i <input file> -c <configuration file> -o <output file> -update <Mean Frechet "
                  << "or Mean Vector> -assignment <Classic or LSH or Hypercube or LSH_Frechet> "
                  << "-complete <optional> -silhouette <optional>" << std::endl;
        return -1;
    }
    params.print_NN_params();

    // read the datasets as vector of Items (Item is described in utils.hpp)
    std::vector<Item> *dataset = new vector<Item>;
    read_items(dataset, params.input_f);

    // convert it to dataset of curves
    vector<curves::Curve2d> *curves_dataset = new vector<curves::Curve2d>;
    // create a vector that will help us represent time
    vector<double> t_dimension;
    for (int i = 0; i < (*dataset)[0].xij.size(); i++)
    {
        t_dimension.push_back(i);
    }

    // create a dataset of curves using our original dataset and the time vector
    for (int i = 0; i < 100; i++)
    {
        curves_dataset->push_back(curves::Curve2d((*dataset)[i].id, t_dimension, (*dataset)[i].xij));
    }

    if (lc(params.update) == "mean vector")
    {
        init::VectorInitializer initializer = init::VectorInitializer(params.clusters, (*dataset), (*dataset)[0].xij.size());
        std::vector<Item> centroids = initializer.initialize_pp();
    }
    if (lc(params.update) == "mean frechet")
    {
        std::cout << "Initializing..." << std::endl;
        init::CurveInitializer initializer = init::CurveInitializer(params.clusters, (*curves_dataset));
        std::vector<curves::Curve2d> centroids = initializer.initialize_pp();
        for (int i = 0; i < centroids.size(); i++)
        {
            std::cout << "id" << i << ": " << centroids[i].id << ", ";
        }
    }

    return 0;
}
