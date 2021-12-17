#ifndef LSH_HPP
#define LSH_HPP
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include "./hashing.hpp"

// class to hold parameters of LSH algorithm
class LSH_params
{
public:
    string input_file = "./datasets/input_small_id";
    string query_file = "./datasets/query_small_id";
    string out_file = "./";
    int k = 4;                 // number of hi functions
    int L = 5;                 // number of hash tables
    int N = 1;                 // number of nearest neighbors to be found
    double R = 10000;          // search range
    string source = "default"; // shows if parameters are default or if they chenged (values: "default", "custom")
};

class G
{
    int k; // number of H functions used in the algorithm
    int tableSize;
    int w; // windowSize inserted by the user
    int d; // dimension of vector
    long unsigned int m;
    default_random_engine eng;
    uniform_int_distribution<int> uid;
    vector<int> r;

public:
    G(int k, int tableSize, int window, int dimensions) : k(k), tableSize(tableSize), w(window), d(dimensions), eng(time(0) + clock()), uid(0, w - 1)
    {
        m = (long unsigned int)(((long long)1 << 32) - (long long)5);
        // create a vector v whose points follow the uniform real distribution
        for (int i = 0; i < k; i++)
        {
            this->r.push_back(uid(eng));
        }
    };

    unsigned int produce_g(const Item &p)
    {
        H h_familly = H(w, d, k);
        vector<int> h = h_familly.produce_k_h(p);
        // for (int j = 0; j < k; j++)
        // {
        //     cout << h[j] << ", ";
        // }
        // cout << endl;

        long unsigned int sum = 0;
        for (int i = 0; i < k; i++)
        {
            /*
            We are using the following formula to not achieve overflow:
            (a$b)modM = ((amodM)$(bmodM))modM.
            $ represents {+, -, *, /}
            */
            // Calculate a mod M
            long unsigned int amodm = (long unsigned)this->r[i] % m;

            // Calculate b mod M
            long unsigned int bmodm = (long unsigned)h[i] % m;

            // Calculate (amodM*bmodM)modM
            long unsigned int resmodm = (amodm * bmodm) % m;

            // Add the above into the sum
            sum += resmodm;
        }
        sum = sum % m;

        // sum = sum % (long unsigned)tableSize;
        return (unsigned int)sum;
    }
};

class LSH
{
    int dimension; // Data dimension

    int windowSize; // windowSize

    int tableSize; // hashtableSize

    LSH_params params; // k, L, N, R

    std::vector<Item *> **hashTables;
    G **g;

public:
    LSH(const LSH_params &params, vector<Item> &dataset, int factor_for_windowSize, int divisor_for_tableSize) : params(params)
    {
        // tune windowSize
        tableSize = dataset.size() / divisor_for_tableSize;
        dimension = dataset[0].xij.size();

        std::random_device rd;                                         // only used once to initialise (seed) engine
        std::mt19937 rng(rd());                                        // random-number engine used (Mersenne-Twister in this case)
        std::uniform_int_distribution<int> uni(0, dataset.size() - 1); // guaranteed unbiased
        int item_index_1;
        int item_index_2;
        double distance = 0;

        /* For dataset.size()/4 samples we randomly choose two points of the dataset and calculate their Euclidean distance.
            We sum these disances and calculate the average. Then we multiply that average distance by a factor of our choice and we
            get the windowsize. That factor can be adjusted at LSH.cpp at the initialization of the LSH instance */
        for (int i = 0; i < dataset.size() / 4; i++)
        {
            item_index_1 = uni(rng);
            item_index_2 = uni(rng);
            while (item_index_1 == item_index_2)
                item_index_2 = uni(rng);
            distance += (EuclideanDistance(&dataset[item_index_1], &dataset[item_index_2], dimension)) / (double)(dataset.size() / 4);
        }

        windowSize = factor_for_windowSize * (int)distance;
        // cout << "w " << w << endl;

        // Initialize L hashTables and g_hashFunctions
        hashTables = new std::vector<Item *> *[params.L];
        g = new G *[params.L];
        for (int i = 0; i < params.L; i++) // for every hashTable
        {
            hashTables[i] = new std::vector<Item *>[tableSize];
            g[i] = new G(params.k, tableSize, windowSize, dimension);
            // cout << g[i]->produce_g(dataset[0]) << " ";
        }
        // cout << endl;

        // Hash all items in training set and insert them into their buckets
        for (int a = 0; a < dataset.size(); a++)
        {
            for (int i = 0; i < params.L; i++)
            {
                unsigned int bucket = g[i]->produce_g(dataset[a]) % (long unsigned)tableSize;
                hashTables[i][bucket].push_back(&dataset[a]);
            }
        }
    }
    ~LSH()
    {
        for (int i = 0; i < params.L; i++) // for every hashTable
        {
            delete[] hashTables[i];
            delete g[i];
        }
        delete[] hashTables;
        delete[] g;
    }

    std::vector<std::pair<double, Item *>> kNN(const Item *query, int thresh = 0) const
    {
        int N = params.N;
        // initialize a vector of N best candidates and distances represented as c++ pairs
        std::vector<std::pair<double, Item *>> knns;
        // Then initialize each pair with max distance and a null item
        for (int i = 0; i < N; i++)
        {
            Item item = Item("null");
            knns.push_back(std::make_pair(std::numeric_limits<double>::max(), &item));
        }

        // For each hash table...
        int itemsSearched = 0;
        for (int i = 0; i < params.L; i++)
        {
            // Calculate the bucket to which the query item corresponds
            long unsigned id = g[i]->produce_g(*query);
            long unsigned bucket = id % (long unsigned)tableSize;

            // For each item inside the bucket...
            for (int j = 0; j < hashTables[i][bucket].size(); j++)
            {
                /*
                Check if the current item is already inserted into the KNNs vector from a previous hash table.
                */
                bool alreadyExists = false;
                for (int a = 0; a < N; a++)
                    if (knns[a].second->id == hashTables[i][bucket][j]->id)
                        alreadyExists = true;

                if (alreadyExists)
                    continue;

                // // // Querying trick (from Lecture Slides): ID is locality sensitive. Avoid computing Euclidean distance for all elements in bucket.
                // if(g[i]->produce_g(*(hashTables[i][bucket][j]))!=id)
                //     continue;

                // Calculate item's distance to the query item
                double distance = EuclideanDistance(query, hashTables[i][bucket][j], dimension);

                /*
                The last pair in the N-sized vector is the worst out of the N
                best candidates till now. If a better candidate is found,
                replace the last pair with the new one and re-sort the vector.
                */
                if (distance < knns[N - 1].first)
                {
                    knns[N - 1].first = distance;
                    knns[N - 1].second = hashTables[i][bucket][j];
                    std::sort(knns.begin(), knns.end(), comparePairs);
                }

                /*
                If a certain threshold of items traversed is reached, return the vector.
                If thresh == 0 it indicates that the user does not want to add a threshold.
                */
                itemsSearched++;
                if (thresh != 0 && itemsSearched >= thresh)
                    return knns;
            }
        }
        return knns;
    }

    /*
    Each neighbor is represented as a pair of <distanceToQuery, neighborItem*>
    The following function returns a vector of these pairs
    */
    std::vector<std::pair<double, Item *>> RangeSearch(const Item *query, double radius, int thresh = 0) const
    {
        std::vector<std::pair<double, Item *>> d;

        // For each hash table...
        int itemsSearched = 0;
        for (int i = 0; i < params.L; i++)
        {
            // Calculate the bucket to which the query item corresponds
            long unsigned id = g[i]->produce_g(*query);
            long unsigned bucket = id % (long unsigned)tableSize;

            // For each item inside the bucket...
            for (int j = 0; j < hashTables[i][bucket].size(); j++)
            {
                // Check if the current item is already inserted into the vector
                bool alreadyExists = false;
                for (int a = 0; a < d.size(); a++)
                    if (d[a].second->id == hashTables[i][bucket][j]->id)
                        alreadyExists = true;

                /* In the "reverse assignment with range search using LSH" clustering algorithm we mark items when they are
                assigned to a cluster so the next range search doesn't check them. In ANN all items are unmarked so this
                has no effect */
                if (alreadyExists || hashTables[i][bucket][j]->marked == true)
                    continue;

                double distance = EuclideanDistance(query, hashTables[i][bucket][j], dimension);

                // If the distance is less than radius, insert the pair into the return vector
                if (distance < radius)
                {
                    std::pair<double, Item *> tmpPair = std::make_pair(distance, hashTables[i][bucket][j]);
                    d.push_back(tmpPair);
                }

                // If a certain threshold of items traversed is reached, return the vector.
                itemsSearched++;
                if (thresh != 0 && itemsSearched >= thresh)
                {
                    // std::sort(d.begin(), d.end(), comparePairs);
                    return d;
                }
            }
        }
        // std::sort(d.begin(), d.end(), comparePairs);
        return d;
    }
};

#endif
