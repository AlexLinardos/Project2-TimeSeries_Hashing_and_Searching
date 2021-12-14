#include <iostream>
#include "NN_interface.hpp"
#include "./Basic/LSH.hpp"
#include "./Basic/HC.hpp"
#include "./ContinuousFrechet/con_Frechet.hpp"
#include "./DiscreteFrechet/disc_Frechet.hpp"

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

    // read the datasets as vector of Items (Item is described in utils.hpp)
    vector<Item>* dataset = new vector<Item>;
    read_items(dataset, params.input_f);
    vector<Item>* queries = new vector<Item>;
    read_items(queries, params.query_f);

    if (params.algorithm == "LSH")
    { // pass parameters to LSH_params class so we can use code from previous project
        LSH_params lsh_params;
        lsh_params.input_file = params.input_f;
        lsh_params.query_file = params.query_f;
        lsh_params.k = params.k;
        lsh_params.L = params.L;
        lsh_params.out_file = params.output_f;

        std::cout << "------[LSH]------" << std::endl;
        LSH *lsh = new LSH(lsh_params, dataset, 3, 8);

        std::cout << "Searching for " << lsh_params.N << " nearest neighbors..." << std::endl;

        std::vector<std::pair<double, Item *>> knns;
        std::vector<std::pair<double, Item *>> true_knns;
        std::vector<std::pair<double, Item *>> r;

        std::chrono::steady_clock::time_point lsh_begin;
        std::chrono::steady_clock::time_point lsh_end;
        std::chrono::steady_clock::time_point true_begin;
        std::chrono::steady_clock::time_point true_end;

        double error = 0.0;
        double avg_error = 0.0;

        ofstream output_file;
        output_file.open(params.output_f);
        double lsh_elapsed = 0.0;
        double brute_elapsed = 0.0;
        clock_t begin;
        clock_t end;

        for (int i = 0; i < queries.size(); i++)
        {
            output_file << "Query: " << queries[i].id << endl;

            // cout << "[k-ANN]" << endl;
            lsh_begin = std::chrono::steady_clock::now();
            begin = clock();
            knns = lsh->kNN(&queries[i], dataset.size() / 5);
            end = clock();
            lsh_end = std::chrono::steady_clock::now();
            lsh_elapsed += double(end - begin);

            // cout << "[Brute Force]" << endl;
            true_begin = std::chrono::steady_clock::now();
            begin = clock();
            true_knns = brute_force_search(dataset, &queries[i], lsh_params.N);
            end = clock();
            true_end = std::chrono::steady_clock::now();
            brute_elapsed += double(end - begin);
            int neighboors_returned = 0;

            for (int j = 0; j < lsh_params.N; j++)
            {
                if (knns[j].second->null == true)
                {
                    output_file << "Nearest neighbor-" << j + 1 << ": "
                                << "NOT FOUND" << endl;
                    continue;
                }
                output_file << "Nearest neighbor-" << j + 1 << ": " << knns[j].second->id << endl;
                output_file << "distanceLSH: " << knns[j].first << endl;
                output_file << "distanceTrue: " << true_knns[j].first << endl;
                error += (knns[j].first / true_knns[j].first); // sum distLSH/distTrue of the N nearest neigbors of a Query
                neighboors_returned++;
            }
            output_file << "tLSH: " << (std::chrono::duration_cast<std::chrono::microseconds>(lsh_end - lsh_begin).count()) / 1000000.0 << std::endl;
            output_file << "tTrue: " << (std::chrono::duration_cast<std::chrono::microseconds>(true_end - true_begin).count()) / 1000000.0 << std::endl;

            output_file << "R-near neighbors:" << endl;
            r = lsh->RangeSearch(&queries[i], lsh_params.R, dataset.size() / 5);
            for (int a = 0; a < r.size(); a++)
            {
                output_file << r[a].second->id << ", " << r[a].first << endl;
            }

            if (error != 0)
            {
                avg_error += error / (double)neighboors_returned; // sum avg(distLSH/distTrue) of the N nearest neigbors of a Query for all the Queries
                error = 0;
            }
            neighboors_returned = 0;
        }

        cout << "[EVALUATION]" << endl;

        lsh_elapsed = lsh_elapsed / CLOCKS_PER_SEC;
        brute_elapsed = brute_elapsed / CLOCKS_PER_SEC;
        cout << "tlSH/tTrue: " << (lsh_elapsed * 1000000.0) / (brute_elapsed * 1000000.0) << endl;

        avg_error = avg_error / (double)queries.size(); // calculate avg distLSH/distTrue across all Queries

        cout << "distLSH/distTrue (avg): " << avg_error << endl;

        output_file.close();

        delete lsh;
    }

    if (params.algorithm == "Hypercube")
    { // pass parameters to Cube_params class so we can use code from previous project
        Cube_params cube_params;
        cube_params.input_file = params.input_f;
        cube_params.query_file = params.query_f;
        cube_params.k = params.k;
        cube_params.M = params.M;
        cube_params.probes = params.probes;
        cube_params.out_file = params.output_f;

        std::cout << "------[Hypercube]------" << std::endl;

        F f = F(params.k);
        Hypercube *cube = new Hypercube(cube_params, dataset, 3, f.h_maps);

        std::cout << "Searching for " << cube_params.N << " nearest neighbors..." << std::endl;

        std::vector<std::pair<double, Item *>> knns;
        std::vector<std::pair<double, Item *>> true_knns;
        std::vector<std::pair<double, Item *>> r;

        std::chrono::steady_clock::time_point lsh_begin;
        std::chrono::steady_clock::time_point lsh_end;
        std::chrono::steady_clock::time_point true_begin;
        std::chrono::steady_clock::time_point true_end;

        double error = 0.0;
        double avg_error = 0.0;

        ofstream output_file;
        output_file.open(cube_params.out_file);
        double lsh_elapsed = 0.0;
        double brute_elapsed = 0.0;
        clock_t begin;
        clock_t end;

        for (int i = 0; i < queries.size(); i++)
        {
            output_file << "Query: " << queries[i].id << endl;

            // cout << "[k-ANN]" << endl;
            lsh_begin = std::chrono::steady_clock::now();
            begin = clock();
            knns = cube->kNN(&queries[i]);
            end = clock();
            lsh_end = std::chrono::steady_clock::now();
            lsh_elapsed += double(end - begin);

            // cout << "[Brute Force]" << endl;
            true_begin = std::chrono::steady_clock::now();
            begin = clock();
            true_knns = brute_force_search(dataset, &queries[i], cube_params.N);
            end = clock();
            true_end = std::chrono::steady_clock::now();
            brute_elapsed += double(end - begin);
            int neighboors_returned = 0;

            for (int j = 0; j < cube_params.N; j++)
            {
                if (knns[j].second->null == true)
                {
                    output_file << "Nearest neighbor-" << j + 1 << ": "
                                << "NOT FOUND" << endl;
                    continue;
                }
                output_file << "Nearest neighbor-" << j + 1 << ": " << knns[j].second->id << endl;
                output_file << "distanceCUBE: " << knns[j].first << endl;
                output_file << "distanceTrue: " << true_knns[j].first << endl;
                error += (knns[j].first / true_knns[j].first); // sum distLSH/distTrue of the N nearest neigbors of a Query
                neighboors_returned++;
            }
            output_file << "tCUBE: " << (std::chrono::duration_cast<std::chrono::microseconds>(lsh_end - lsh_begin).count()) / 1000000.0 << std::endl;
            output_file << "tTrue: " << (std::chrono::duration_cast<std::chrono::microseconds>(true_end - true_begin).count()) / 1000000.0 << std::endl;

            output_file << "R-near neighbors:" << endl;
            r = cube->RangeSearch(&queries[i], cube_params.R);
            for (int a = 0; a < r.size(); a++)
            {
                output_file << r[a].second->id << ", " << r[a].first << endl;
            }
            if (error != 0)
            {
                avg_error += error / (double)neighboors_returned; // sum avg(distLSH/distTrue) of the returned nearest neigbors of a Query for all the Queries
                error = 0;
            }
            neighboors_returned = 0;
            // cout << meso_error << endl;
        }

        cout << "[EVALUATION]" << endl;

        lsh_elapsed = lsh_elapsed / CLOCKS_PER_SEC;
        brute_elapsed = brute_elapsed / CLOCKS_PER_SEC;
        cout << "tCUBE/tTrue: " << (lsh_elapsed * 1000000.0) / (brute_elapsed * 1000000.0) << endl;

        avg_error = avg_error / (double)queries.size(); // calculate avg distLSH/distTrue across all Queries

        cout << "distCUBE/distTrue (avg): " << avg_error << endl;

        output_file.close();

        delete cube;
    }

    if (params.algorithm == "Frechet")
    {
        std::cout << "------[" << params.metric << " Frechet]------" << std::endl;
        std::cout << "RESULT: " << dF::discrete_frechet(dataset[0], dataset[1]);
    }

    return 0;
}
