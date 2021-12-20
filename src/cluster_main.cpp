#include <iostream>
#include <string>
#include <vector>
#include "../ui/Clustering_interface.hpp"
#include "../TimeSeries-Clustering/initialization.hpp"
#include "../TimeSeries-Clustering/assignment.hpp"
#include "../TimeSeries-Clustering/update.hpp"
#include "../TimeSeries-ANN/L2/LSH.hpp"
#include "../TimeSeries-ANN/L2/HC.hpp"

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

    if (lc(params.update) == "mean vector")
    {
        // INITIALIZATION STEP
        std::cout << "Initializing centroids..." << std::endl;
        init::VectorInitializer initializer = init::VectorInitializer(params.clusters, (*dataset), (*dataset)[0].xij.size());
        std::vector<Item> centroids = initializer.initialize_pp();

        // ASSIGNMENT STEP
        std::cout << "Clustering..." << std::endl;
        assign::VectorAssignor assignor = assign::VectorAssignor(params.clusters, centroids, (*dataset), (*dataset)[0].xij.size());
        if (lc(params.assignment) == "classic")
        {
            assignor.Lloyds();
            // print results
            for (int i = 0; i < assignor.clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < assignor.clusters[i].size(); j++)
                {
                    std::cout << assignor.clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }
        }
        else if (lc(params.assignment) == "lsh")
        {
            // we must construct an LSH object to pass as parameter in the assignment algorithm
            LSH_params lsh_params;
            lsh_params.k = params.k;
            lsh_params.L = params.L;

            LSH lsh_object = LSH(lsh_params, (*dataset), 1.0, 8);
            assignor.Range_LSH(lsh_object);
            // print results
            for (int i = 0; i < assignor.clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < assignor.clusters[i].size(); j++)
                {
                    std::cout << assignor.clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }
        }
        else if (lc(params.assignment) == "hypercube")
        {
            // we must construct a Hypercube object to pass as parameter in the assignment algorithm
            Cube_params hc_params;
            hc_params.k = params.k;
            hc_params.M = params.M;
            hc_params.probes = params.probes;

            F f = F(hc_params.k);
            Hypercube hc_object = Hypercube(hc_params, (*dataset), 1.0, f.h_maps);
            assignor.Range_HC(hc_object);
            // print results
            for (int i = 0; i < assignor.clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < assignor.clusters[i].size(); j++)
                {
                    std::cout << assignor.clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }
        }

        // UPDATE STEP
        std::cout << "------------------------------------" << endl;
        std::cout << "Current centers are: " << endl;
        for (int i = 0; i < assignor.centers.size(); i++)
        {
            for (int j = 0; j < 10; j++)
            {
                std::cout << assignor.centers[i].xij[j] << " ";
            }
            std::cout << std::endl;
        }
        updt::update_vector_centers((*dataset), assignor.centers, assignor.clusters);
        std::cout << "------------------------------------" << endl;
        std::cout << "Updated centers. New centers are: " << endl;
        // print results
        for (int i = 0; i < assignor.centers.size(); i++)
        {
            for (int j = 0; j < 10; j++)
            {
                std::cout << assignor.centers[i].xij[j] << " ";
            }
            std::cout << std::endl;
        }
    }
    if (lc(params.update) == "mean frechet")
    {
        // convert dataset to curves
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

        // INITIALIZATION STEP
        std::cout << "Initializing centroids..." << std::endl;
        init::CurveInitializer initializer = init::CurveInitializer(params.clusters, (*curves_dataset));
        std::vector<curves::Curve2d> centroids = initializer.initialize_pp();

        // ASSIGNMENT STEP
        std::cout << "Clustering..." << std::endl;
        assign::CurveAssignor assignor = assign::CurveAssignor(params.clusters, centroids, (*curves_dataset));
        if (lc(params.assignment) == "classic")
        {
            assignor.Lloyds();
            // print results
            for (int i = 0; i < assignor.clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < assignor.clusters[i].size(); j++)
                {
                    std::cout << assignor.clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }
        }
        else if (lc(params.assignment) == "lsh_frechet")
        {
            double delta = delta_tuning(*curves_dataset);

            // we must construct a discrete Frechet LSH object to pass as parameter in the assignment algorithm
            dFLSH::LSH *dflsh_object = new dFLSH::LSH(curves_dataset, params.L, delta, 8);

            assignor.Range_dfLSH(*dflsh_object);
            // print results
            for (int i = 0; i < assignor.clusters.size(); i++)
            {
                std::cout << "CLUSTER " << i << std::endl;
                for (int j = 0; j < assignor.clusters[i].size(); j++)
                {
                    std::cout << assignor.clusters[i][j].id << " ";
                }
                std::cout << std::endl;
            }
            delete dflsh_object;
        }

        // UPDATE STEP
    }

    return 0;
}
