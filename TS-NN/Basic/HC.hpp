#ifndef HC_HPP
#define HC_HPP
#include <iostream>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <typeinfo>
#include <chrono>
#include <vector>
#include "hashing.hpp"
#include "utils.hpp"

using namespace std;

// class to hold parameters of Hypercube algorithm
class Cube_params
{
public:
    string input_file = "./datasets/input_small_id";
    string query_file = "./datasets/query_small_id";
    string out_file = "./";
    int k = 14;                // dimension at which points will be projected
    int M = 10;                // maximum number of points to be checked
    int probes = 2;            // maximum number of vertices(probes) to be checked
    int N = 1;                 // number of nearest neighbors to be found
    double R = 1000;           // search range
    string source = "default"; // shows if parameters are default or if they chenged (values: "default", "custom")
};

// function that assigns values produced from H family functions to {0,1} uniformly
class F
{
public:
    int k; // dimension at which points will be projected
    vector<unordered_map<int, int>> h_maps;
    default_random_engine eng;
    uniform_int_distribution<int> uid;

    F(int k) : k(k), h_maps(k), eng(chrono::system_clock::now().time_since_epoch().count()), uid(0, 1){};

    int produce_f(int index, int h)
    {
        if (h_maps[index].count(h) == 0)
        {
            int rnum = uid(eng);
            h_maps[index][h] = rnum;
            // cout << "h" << index << ":" << h << "-> " << h_maps[index][h] << endl;
            return rnum;
        }
        else
        {
            return h_maps[index].find(h)->second;
        }
    }
};

class Hypercube
{
public:
    int k;                 // dimension at which points will be projected
    int M;                 // maximum number of points to be checked
    int probes;            // maximum number of vertices(probes) to be checked
    int N;                 // number of nearest neighbors to be found
    double R;              // search range
    int d;                 // dimensions of points (items)
    int w;                 // Window
    unsigned int vertices; // Hypercube vertices

    std::vector<Item *> *hash_table;
    F f;

    Hypercube(const Cube_params &params, vector<Item> &dataset, int factor_for_windowSize, vector<unordered_map<int, int>> &h_maps) : f(params.k)
    {
        d = dataset[0].xij.size();
        std::random_device rd;                                         // only used once to initialise (seed) engine
        std::mt19937 rng(rd());                                        // random-number engine used (Mersenne-Twister in this case)
        std::uniform_int_distribution<int> uni(0, dataset.size() - 1); // guaranteed unbiased
        int item_index_1;
        int item_index_2;
        double distance = 0;
        double length = (double)(dataset.size() / 4);

        // calculate window w according to our dataset
        for (int i = 0; i < length; i++)
        {
            item_index_1 = uni(rng);
            item_index_2 = uni(rng);
            while (item_index_1 == item_index_2)
                item_index_2 = uni(rng);
            distance += EuclideanDistance(&dataset[item_index_1], &dataset[item_index_2], d) / length;
        }

        w = factor_for_windowSize * distance;

        this->k = params.k;
        this->M = params.M;
        this->probes = params.probes;
        this->N = params.N;
        this->R = params.R;
        this->vertices = pow(2, k);

        hash_table = new std::vector<Item *>[vertices];

        // hash_table = vector<vector<Item>>(vertices);
        // cout << "INITIALIZED HASH TABLE OF SIZE: " << hash_table.size() << endl;
        // for (int i = 0; i < vertices; ++i)
        // {
        //     vector<Item> dummy_vec;
        //     hash_table[i] = dummy_vec;
        // }

        for (int i = 0; i < dataset.size(); ++i)
        {
            // produce the needed H hash functions for each point
            H h(w, d, k);
            vector<int> h_functions = h.produce_k_h(dataset[i]);
            vector<int> f_values(k);

            // assign H values to {0,1} using F function
            for (int j = 0; j < k; ++j)
            {
                f_values.push_back(f.produce_f(j, h_functions[j]));
            }
            // concatenate f values to create a bitstring
            unsigned int bucket = concat_f_values(f_values);
            // assign item to the corresponding bucket using the above bitstring
            hash_table[bucket].push_back(&dataset[i]);
        }
    }

    ~Hypercube()
    {
        delete[] hash_table;
    }

    // calculates the bucket (vertex) at which query should be assigned to
    unsigned int find_bucket(const Item &q)
    {
        H h(w, d, k);
        vector<int> h_functions = h.produce_k_h(q);
        vector<int> f_values(k);
        for (int j = 0; j < k; ++j)
        {
            f_values.push_back(f.produce_f(j, h_functions[j]));
        }
        return concat_f_values(f_values);
    }

    // concatenates a series of values {0,1} using them as bits to create an unsigned integer
    unsigned int concat_f_values(const vector<int> &f_values) const
    {
        unsigned int num = 0;
        for (int i = 0; i < f_values.size(); ++i)
        {
            if (f_values[i] == 1)
            {
                num |= 1;
            }

            if (i < (f_values.size() - 1))
            {
                num <<= 1;
            }
        }
        return num;
    }

    // returns probes at Hamming distance dist of given bucket (vertex)
    vector<int> get_probes_in_distance(int dist, unsigned int bucket) const
    {
        vector<int> result;

        for (int i = 0; i < this->vertices; ++i)
        {
            if (HammingDistance(bucket, i) == dist)
            {
                result.push_back(i);
            }
        }
        return result;
    }

    // returns all probes until it reaches threshold parameter
    // begins from the closests and continues by incrementing the distance if needed
    vector<int> get_probes_in_threshold(unsigned int bucket) const
    {
        vector<int> result;

        int probes_reached = 0;
        for (int i = 0; i < this->vertices; ++i)
        {
            vector<int> p_in_dist = get_probes_in_distance(i, bucket);
            // cout << "DISTANCE" << i << ": ";
            for (int j = 0; j < p_in_dist.size(); ++j)
            {
                // cout << p_in_dist[j] << " ";
            }
            // cout << endl;

            for (int j = 0; j < p_in_dist.size(); ++j)
            {
                result.push_back(p_in_dist[j]);
                probes_reached++;
                if (probes_reached >= this->probes)
                {
                    // cout << "Ended early at " << probes_reached << " probes." << endl;
                    return result;
                }
            }
        }
        // cout << "I reached " << probes_reached << " instead of " << this->probes << endl;
        return result;
    }

    // performs the kNN algorithm using Hypercube randomized projection
    std::vector<std::pair<double, Item *>> kNN(Item *query)
    {
        // At first initalize the result vector of <distanceFromQuery, item> pairs
        std::vector<std::pair<double, Item *>> knns;
        // Then initialize each pair with distance -> (max integer) and a null item

        for (int i = 0; i < N; i++)
        {
            Item item = Item("null");
            knns.push_back(std::make_pair(std::numeric_limits<double>::max(), &item));
        }

        unsigned int q_bucket = find_bucket(*query);
        // cout << "Query " << query.id << " is in bucket " << q_bucket << endl;
        vector<int> search_probes = get_probes_in_threshold(q_bucket);

        // cout << "-----------------------------------------------------" << endl;
        int items_searched = 0;
        for (int i = 0; i < search_probes.size(); ++i)
        {
            int curr_bucket = search_probes[i];
            // cout << "Searching bucket " << curr_bucket << " with " << hash_table[curr_bucket].size() << " items inside." << endl;
            // cout << "BUCKET " << curr_bucket << ": ";
            //  for (int j = 0; j < hash_table[curr_bucket].size(); ++j)
            //  {
            //      cout << hash_table[curr_bucket][j].id << ", ";
            //  }
            //  cout << endl;

            for (int j = 0; j < hash_table[curr_bucket].size(); ++j)
            {
                // cout << "Item " << hash_table[curr_bucket][j].id << " in ";
                double dist = EuclideanDistance(query, hash_table[curr_bucket][j], d);
                // cout << "distance " << dist << " | ";

                if (dist < knns[N - 1].first)
                {
                    knns[N - 1].first = dist;
                    knns[N - 1].second = hash_table[curr_bucket][j];
                    std::sort(knns.begin(), knns.end(), comparePairs);
                }

                items_searched++;
                if (items_searched >= this->M)
                {
                    return knns;
                }
            }
            // cout << endl;
            // cout << "...................." << endl;
        }
        return knns;
    }

    // performs the Range Search algorithm using Hypercube randomized projection
    std::vector<std::pair<double, Item *>> RangeSearch(Item *query, double radius)
    {
        // At first initalize the result vector of <distanceFromQuery, item> pairs
        std::vector<std::pair<double, Item *>> rns;
        unsigned int q_bucket = find_bucket(*query);
        // cout << "Query " << query.id << " is in bucket " << q_bucket << endl;
        vector<int> search_probes = get_probes_in_threshold(q_bucket);

        // cout << "-----------------------------------------------------" << endl;
        int items_searched = 0;
        int curr_bucket;
        for (int i = 0; i < search_probes.size(); ++i)
        {
            curr_bucket = search_probes[i];
            // cout << "Searching bucket " << curr_bucket << " with " << hash_table[curr_bucket].size() << " items inside." << endl;
            // cout << "BUCKET " << curr_bucket << ": ";
            // for (int j = 0; j < hash_table[curr_bucket].size(); ++j)
            // {
            //     cout << hash_table[curr_bucket][j].id << ", ";
            // }
            // cout << endl;

            for (int j = 0; j < hash_table[curr_bucket].size(); ++j)
            {
                /* In the "reverse assignment with range search" clustering algorithm we mark items when they are
                assigned to a cluster so the next range search doesn't check them. In ANN all items are unmarked so this
                has no effect */
                if (hash_table[curr_bucket][j]->marked == true)
                    continue;

                // cout << "Item " << hash_table[curr_bucket][j].id << " in ";
                double dist = EuclideanDistance(query, hash_table[curr_bucket][j], d);
                // cout << "distance " << dist << " | ";

                if (dist < radius)
                {
                    std::pair<double, Item *> tmp_pair = std::make_pair(dist, hash_table[curr_bucket][j]);
                    rns.push_back(tmp_pair);
                }

                items_searched++;
                if (items_searched >= this->M)
                {
                    return rns;
                }
            }
            // cout << endl;
            // cout << "...................." << endl;
        }
        return rns;
    }
};

#endif
