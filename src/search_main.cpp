#include <iostream>
#include <vector>
#include <string>
#include "../ui/NN_interface.hpp"
#include "../includes/curves.hpp"
#include "../includes/mean_curve.hpp"
#include "../TimeSeries-ANN/L2/LSH.hpp"
#include "../TimeSeries-ANN/L2/HC.hpp"
#include "../TimeSeries-ANN/DiscreteFrechet/discF_LSH.hpp"
#include "../TimeSeries-ANN/ContinuousFrechet/contF_LSH.hpp"

// ./bin/search -i <input file> -q <query file> -k 1 -L 2 -M 3 -probes 4 -o <output file> -algorithm <LSH> -metric <discrete> -delta 5.01
int main(int argc, char *argv[])
{
    NNi::NN_params params = NNi::NN_params(argc, argv);
    if (params.success == false)
    {
        std::cout << "-----------------------------------------------------------" << std::endl
                  << "Something went wrong while reading command line parameters." << std::endl
                  << "Please make sure you follow the format bellow: " << std::endl
                  << "./bin/search -i <input file> -q <query file> -k <int> -L <int> -M <int> -probes "
                  << "<int> -o <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete "
                  << "or continuous | only for -algorithm Frechet> -delta <double>" << std::endl;
        return -1;
    }
    params.print_NN_params();

    // read the datasets as vector of Items (Item is described in utils.hpp)
    vector<Item> *dataset = new vector<Item>;
    read_items(dataset, params.input_f);
    vector<Item> *queries = new vector<Item>;
    read_items(queries, params.query_f);

    std::chrono::steady_clock::time_point lsh_begin;
    std::chrono::steady_clock::time_point true_begin;

    double error = 0.0;
    double lsh_elapsed = 0.0;
    double brute_elapsed = 0.0;
    double maf = 0.0;
    double f = 0.0;

    if (lc(params.algorithm) == "lsh") // ideal is L=1, k=1, tablesize=datasetsize/2, w = average L2 distance between dataset curves (as vectors)
    { // pass parameters to Cube_params class so we can use code from previous project
        LSH_params lsh_params;
        lsh_params.input_file = params.input_f;
        lsh_params.query_file = params.query_f;
        lsh_params.k = params.k;
        lsh_params.L = params.L;
        lsh_params.out_file = params.output_f;
        lsh_params.N = 1;

        // lsh_params.k = 1;
        // lsh_params.L = 1; 

        std::vector<std::pair<double, Item *>> knns;
        std::vector<std::pair<double, Item *>> true_knns;

        std::cout << "------[LSH]------" << std::endl;

        LSH *lsh = new LSH(lsh_params, *dataset, 1, 2);

        std::cout << "Searching for the approximate nearest neighbors of the query curves..." << std::endl;

        ofstream output_file;
        output_file.open(params.output_f);

        output_file << "Algorithm:  LSH_Vector" << endl
                    << endl;
        for (int i = 0; i < queries->size(); i++)
        {
            output_file << "Query: " << (*queries)[i].id << endl;

            // cout << "[ANN]" << endl;
            lsh_begin = std::chrono::steady_clock::now();
            knns = lsh->kNN(&(*queries)[i], dataset->size() / 4);
            lsh_elapsed += (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - lsh_begin).count();

            // cout << "[Brute Force]" << endl;
            true_begin = std::chrono::steady_clock::now();
            true_knns = brute_force_search((*dataset), &((*queries)[i]), lsh_params.N);
            brute_elapsed += (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - true_begin).count();

            int neighboors_returned = 0;

            if (knns[0].second->null == true)
            {
                output_file << "Approximate Nearest neighbor "
                            << "NOT FOUND" << endl;
                continue;
            }
            output_file << "Approximate Nearest neighbor: " << knns[0].second->id << endl;
            output_file << "True Nearest neighbor: " << true_knns[0].second->id << endl;
            output_file << "distanceApproximate: " << knns[0].first << endl;
            output_file << "distanceTrue: " << true_knns[0].first << endl;
            error += (knns[0].first / true_knns[0].first); // sum distLSH/distTrue of the nearest neigbor of a Query

            neighboors_returned++;

            // υπολογίζουμε dist(approx NN) / dist(true NN) για κάθε query και κρατάμε το max όλων
            f = knns[0].first / true_knns[0].first;
            if (f > maf)
                maf = f;
            output_file << endl;
        }
        output_file << endl
                    << "tApproximateAverage: " << lsh_elapsed / (double)queries->size() << " (μs)" << endl;
        output_file << "tTrueAverage: " << brute_elapsed / (double)queries->size() << " (μs)" << endl;
        output_file << "MAF: " << maf << endl;

        cout << "[EVALUATION]" << endl;
        cout << "tlSH/tTrue: " << lsh_elapsed / brute_elapsed << endl;
        cout << "distLSH/distTrue (avg): " << error / (double)queries->size() << endl;

        output_file.close();

        delete lsh;
    }

    if (lc(params.algorithm) == "hypercube") // ideal is k=1. M=30, probes=1, w = average L2 distance between dataset curves (as vectors)
    { // pass parameters to Cube_params class so we can use code from previous project
        Cube_params cube_params;
        cube_params.input_file = params.input_f;
        cube_params.query_file = params.query_f;
        cube_params.k = params.k;
        cube_params.M = params.M;
        cube_params.probes = params.probes;
        cube_params.out_file = params.output_f;
        cube_params.N = 1;

        cube_params.k = 1;
        cube_params.M = 30;
        cube_params.probes = 1;

        std::vector<std::pair<double, Item *>> knns;
        std::vector<std::pair<double, Item *>> true_knns;

        std::cout << "------[Hypercube]------" << std::endl;

        F f_map = F(params.k);
        Hypercube *cube = new Hypercube(cube_params, *dataset, 1, f_map.h_maps);

        std::cout << "Searching for the approximate nearest neighbors of the query curves..." << std::endl;

        ofstream output_file;
        output_file.open(params.output_f);

        output_file << "Algorithm:  Hypercube" << endl
                    << endl;
        for (int i = 0; i < queries->size(); i++)
        {
            output_file << "Query: " << (*queries)[i].id << endl;

            // cout << "[ANN]" << endl;
            lsh_begin = std::chrono::steady_clock::now();
            knns = cube->kNN(&(*queries)[i]);
            lsh_elapsed += (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - lsh_begin).count();

            // cout << "[Brute Force]" << endl;
            true_begin = std::chrono::steady_clock::now();
            true_knns = brute_force_search(*dataset, &(*queries)[i], cube_params.N);
            brute_elapsed += (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - true_begin).count();

            int neighboors_returned = 0;

            if (knns[0].second->null == true)
            {
                output_file << "Approximate Nearest neighbor "
                            << "NOT FOUND" << endl;
                continue;
            }
            output_file << "Approximate Nearest neighbor: " << knns[0].second->id << endl;
            output_file << "True Nearest neighbor: " << true_knns[0].second->id << endl;
            output_file << "distanceApproximate: " << knns[0].first << endl;
            output_file << "distanceTrue: " << true_knns[0].first << endl;
            error += (knns[0].first / true_knns[0].first); // sum distLSH/distTrue of the nearest neigbor of a Query

            neighboors_returned++;

            // υπολογίζουμε dist(approx NN) / dist(true NN) για κάθε query και κρατάμε το max όλων
            f = knns[0].first / true_knns[0].first;
            if (f > maf)
                maf = f;
            output_file << endl;
        }
        output_file << endl
                    << "tApproximateAverage: " << lsh_elapsed / (double)queries->size() << " (μs)" << endl;
        ;
        output_file << "tTrueAverage: " << brute_elapsed / (double)queries->size() << " (μs)" << endl;
        ;
        output_file << "MAF: " << maf << endl;

        cout << "[EVALUATION]" << endl;
        cout << "tCUBE/tTrue: " << lsh_elapsed / brute_elapsed << endl;
        cout << "distCUBE/distTrue (avg): " << error / (double)queries->size() << endl;

        output_file.close();

        delete cube;
    }

    if (lc(params.algorithm) == "frechet")
    {
        std::cout << "------[" << params.metric << " Frechet]------" << std::endl;

        vector<curves::Curve2d> *curves_dataset = new vector<curves::Curve2d>;
        vector<curves::Curve2d> *filtered_curves_dataset = new vector<curves::Curve2d>;
        vector<curves::Curve2d> *curves_queryset = new vector<curves::Curve2d>;
        vector<curves::Curve2d> *filtered_curves_queryset = new vector<curves::Curve2d>;

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

        // do the same for query set
        for (int i = 0; i < (*queries).size(); i++)
        {
            curves_queryset->push_back(curves::Curve2d((*queries)[i].id, t_dimension, (*queries)[i].xij));
        }

        if(params.delta == 0.0) // if delta is not defined by user or if it is 0.0, auto tune it according to dataset
        {
            params.delta = delta_tuning(*curves_dataset);
        }

        if (lc(params.metric) == "discrete") // L = 6, delta = avg dist between curve vertices, tablesize = dataset_size/8, threshold = dataset_size/4, querying trick = false
        {
            // perform LSH for discrete Frechet

            dFLSH::LSH *dLSH = new dFLSH::LSH(curves_dataset, params.L, params.delta, 8);

            std::cout << "Searching for the approximate nearest neighbors of the query curves..." << std::endl;
            ofstream output_file;
            output_file.open(params.output_f);

            output_file << "Algorithm: LSH_Frechet_Discrete" << endl
                        << endl;
            for (int i = 0; i < curves_queryset->size(); i++)
            {
                output_file << "Query: " << (*curves_queryset)[i].id << endl;

                // cout << "[ANN]" << endl;
                lsh_begin = std::chrono::steady_clock::now();
                std::pair<curves::Curve2d *, double> ann = dLSH->search_ANN((*curves_queryset)[i], false, curves_dataset->size() / 4);
                lsh_elapsed += (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - lsh_begin).count() / 1000000.0);

                // cout << "[Brute Force]" << endl;
                true_begin = std::chrono::steady_clock::now();
                std::pair<curves::Curve2d *, double> true_nn = dF::search_exactNN((*curves_queryset)[i], *curves_dataset);
                brute_elapsed += (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - true_begin).count() / 1000000.0);
                int neighboors_returned = 0;

                if (ann.first->id == "null")
                {
                    output_file << "Approximate Nearest neighbor "
                                << "NOT FOUND" << endl;
                    continue;
                }
                output_file << "Approximate Nearest neighbor: " << ann.first->id << endl;
                output_file << "True Nearest neighbor: " << true_nn.first->id << endl;
                output_file << "distanceApproximate: " << ann.second << endl;
                output_file << "distanceTrue: " << true_nn.second << endl;
                error += (ann.second / true_nn.second); // sum distLSH/distTrue of the nearest neigbor of a Query

                neighboors_returned++;

                // υπολογίζουμε dist(approx NN) / dist(true NN) για κάθε query και κρατάμε το max όλων
                f = ann.second / true_nn.second;
                if (f > maf)
                    maf = f;
                output_file << endl;
            }
            output_file << endl
                        << "tApproximateAverage: " << lsh_elapsed / (double)queries->size() << " (sec)" << endl;
            output_file << "tTrueAverage: " << brute_elapsed / (double)queries->size() << " (sec)" << endl;
            output_file << "MAF: " << maf << endl;

            cout << "[EVALUATION]" << endl;
            cout << "tdF_LSH/tTrue: " << lsh_elapsed / brute_elapsed << endl;
            cout << "distdF_LSH/distTrue (avg): " << error / (double)curves_queryset->size() << endl;

            output_file.close();

            delete dLSH;
        }
        else if (lc(params.metric) == "continuous")
        {
            filtered_curves_dataset = cF::filter_curves(*curves_dataset, 2*params.delta);
            filtered_curves_queryset = cF::filter_curves(*curves_queryset, 2*params.delta);

            // perform LSH for continuous Frechet
            cFLSH::LSH *cLSH = new cFLSH::LSH(filtered_curves_dataset, 1, params.delta, 8);

            std::cout << "Searching for the approximate nearest neighbors of the query curves..." << std::endl;
            ofstream output_file;
            output_file.open(params.output_f);

            output_file << "Algorithm: LSH_Frechet_Continuous" << endl
                        << endl;
            for (int i = 0; i < filtered_curves_queryset->size(); i++)
            {
                output_file << "Query: " << (*curves_queryset)[i].id << endl;

                // cout << "[ANN]" << endl;
                lsh_begin = std::chrono::steady_clock::now();
                std::pair<curves::Curve2d *, double> ann = cLSH->search_ANN((*filtered_curves_queryset)[i], true);
                lsh_elapsed += (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - lsh_begin).count() / 1000000.0);

                // cout << "[Brute Force]" << endl;
                true_begin = std::chrono::steady_clock::now();
                std::pair<curves::Curve2d *, double> true_nn = cF::search_exactNN((*filtered_curves_queryset)[i], *filtered_curves_dataset);
                brute_elapsed += (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - true_begin).count() / 1000000.0);

                int neighboors_returned = 0;

                if (ann.first->id == "null")
                {
                    output_file << "Approximate Nearest neighbor "
                                << "NOT FOUND" << endl;
                    continue;
                }
                output_file << "Approximate Nearest neighbor: " << ann.first->id << endl;
                output_file << "True Nearest neighbor: " << true_nn.first->id << endl;
                output_file << "distanceApproximate: " << ann.second << endl;
                output_file << "distanceTrue: " << true_nn.second << endl;
                error += (ann.second / true_nn.second); // sum distLSH/distTrue of the nearest neigbor of a Query

                neighboors_returned++;

                // υπολογίζουμε dist(approx NN) / dist(true NN) για κάθε query και κρατάμε το max όλων
                f = ann.second / true_nn.second;
                if(f > maf)
                    maf = f;
                output_file << endl; 
            }
            output_file << endl << "tApproximateAverage: " << lsh_elapsed / (double)queries->size() << " (sec)" << endl;
            output_file << "tTrueAverage: " << brute_elapsed / (double)queries->size() << " (sec)" << endl;
            output_file << "MAF: " << maf << endl;

            cout << "[EVALUATION]" << endl;
            cout << "tcF_LSH/tTrue: " << lsh_elapsed / brute_elapsed << endl;
            cout << "distcF_LSH/distTrue (avg): " << error / (double)curves_queryset->size() << endl;

            output_file.close();

            delete cLSH;
        }

        delete filtered_curves_queryset;
        delete filtered_curves_dataset;
        delete curves_dataset;
        delete curves_queryset;
    }

    delete dataset;
    delete queries;
    return 0;
}
