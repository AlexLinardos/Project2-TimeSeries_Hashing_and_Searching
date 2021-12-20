#include <iostream>
#include <string>
#include <vector>
#include "../ui/Clustering_interface.hpp"
#include "../TimeSeries-Clustering/vector_clustering.hpp"
#include "../TimeSeries-Clustering/curve_clustering.hpp"

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

    clock_t begin;
    clock_t end;

    if (lc(params.update) == "mean vector")
    {
        vector_cluster::Clustering *c = new vector_cluster::Clustering(params, dataset);
        if(lc(params.assignment) == "classic")
        {
            begin = clock();
            c->Classic_Clustering(20);
            end = clock();
        }
        else if(lc(params.assignment) == "lsh")
        {
            begin = clock();
            c->Reverse_Assignment_LSH_Clustering(20);
            end = clock();
        }
        else if(lc(params.assignment) == "hypercube")
        {
            begin = clock();
            c->Reverse_Assignment_HC_Clustering(20);
            end = clock();
        }    

        double elapsed = double(end - begin) / CLOCKS_PER_SEC;
        ofstream output_file;
        output_file.open(params.output_f);

        cout << "[BUILDING OUTPUT]" << endl;

        output_file << "Algorithm: " << params.assignment << " " << params.update << endl;

    if (params.complete)
    {
        for (int i = 0; i < params.clusters; i++)
        {
            output_file << "CLUSTER-" << i + 1 << " {centroid: [";
            for (int j = 0; j < (*dataset)[0].xij.size(); j++)
            {
                output_file << c->centers[i].xij[j] << ",";
            }
            output_file << "]";
            for (int j = 0; j < c->clusters[i].size(); j++)
            {
                output_file << ", " << c->clusters[i][j].id;
            }
            output_file << "}" << endl;
        }
    }
    else
    {
        for (int i = 0; i < params.clusters; i++)
        {
            output_file << "CLUSTER-" << i + 1 << " {size: " << c->clusters[i].size() << ", centroid: [";
            for (int j = 0; j < (*dataset)[0].xij.size(); j++)
            {
                output_file << c->centers[i].xij[j] << ",";
            }
            output_file << "]}" << endl;
        }
    }

    if (params.silhouette == true)
    {

        cout << "[CALCULATING SILHOUETTE]" << endl;
        if (c->centers.size() > 1)
        {
            output_file << "clustering_time: " << elapsed << endl;
            output_file << "Silhouette: [";
            for (int i = 0; i < params.clusters; i++)
            {
                // cout << i << endl;
                output_file << c->eval_specific_cluster(i) << ", ";
            }
            output_file << c->eval_clustering() << "]" << endl;
        }
        else
        {
            cout << "Silhouette can only be computed when there are 2 or more clusters." << endl;
        }
    }

        output_file.close();

        delete c;

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

        curve_cluster::Clustering *c = new curve_cluster::Clustering(params, curves_dataset);

        if(lc(params.assignment) == "classic")
        {
            begin = clock();
            c->Classic_Clustering(20);
            end = clock();
        }
        else if(lc(params.assignment) == "lsh_frechet")
        {
            begin = clock();
            c->Reverse_Assignment_LSH_Clustering(20);
            end = clock();
        }

        double elapsed = double(end - begin) / CLOCKS_PER_SEC;
        ofstream output_file;
        output_file.open(params.output_f);

        cout << "[BUILDING OUTPUT]" << endl;

        output_file << "Algorithm: " << params.assignment << " " << params.update << endl;

        if (params.complete == true)
        {
            for (int i = 0; i < params.clusters; i++)
            {
                output_file << "CLUSTER-" << i + 1 << " {centroid: [";
                output_file << "(" << c->centers[i][0].x << "," << c->centers[i][0].y << ")";
                for (int j = 1; j <  c->centers[i].size(); j++)
                {
                    output_file << ", (" << c->centers[i][j].x << "," << c->centers[i][j].y << ")";
                }
                output_file << "]";
                for (int j = 0; j < c->clusters[i].size(); j++)
                {
                    output_file << ", " << c->clusters[i][j].id;
                }
                output_file << "}" << endl;
            }
        }
        else
        {
            for (int i = 0; i < params.clusters; i++)
            {
                output_file << "CLUSTER-" << i + 1 << " {size: " << c->clusters[i].size() << ", centroid: [";
                output_file << "(" << c->centers[i][0].x << "," << c->centers[i][0].y << ")";
                for (int j = 1; j <  c->centers[i].size(); j++)
                {
                    output_file << ", (" << c->centers[i][j].x << "," << c->centers[i][j].y << ")";
                }
                output_file << "]}" << endl;
            }
        }

        if (params.silhouette == true)
        {
            cout << "[CALCULATING SILHOUETTE] (This step takes a couple of minutes)" << endl;
            if (params.clusters > 1)
            {
                output_file << "clustering_time: " << elapsed << endl;
                output_file << "Silhouette: [";
                for (int i = 0; i < params.clusters; i++)
                {
                    // cout << i << endl;
                    output_file << c->eval_specific_cluster(i) << ", ";
                }
                output_file << c->eval_clustering() << "]" << endl;
            }
            else
            {
                cout << "Silhouette can only be computed when there are 2 or more clusters." << endl;
            }
        }

        output_file.close();

        delete c;
    }

    return 0;
}
