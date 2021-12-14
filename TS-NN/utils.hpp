#ifndef UTILITIES_HPP
#define UTILITIES_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <cmath>
#include <unordered_map>
#include <algorithm>

using namespace std;

// each line of our dataset/query file corresponds to an Item with an id and a vector of numeric data
class Item
{
public:
    string id;
    vector<double> xij;

    int cluster = 0;
    int second = -1;

    bool null = false;    // flag if is NULL item created to initialize apair vector for knn and brute foce algorithms
    bool claimed = false; // will be used in reverse assignment to indicate if item has been claimed by a ball (in order to resolve conflicts)
    bool marked = false;  // will be used in reverse assignment to indicate item has been assigned to a cluster

    Item(string id, const vector<double> &xij) : id(id), xij(xij) { null = false; }

    Item(const string &flag)
    {
        if (flag == "null")
        {
            id = "-1";
            null = true;
        }
    }

    Item() {}
};

// used to tokenize a text line
vector<string> tokenize(const string &str)
{
    vector<string> tokens;
    string tmp;
    stringstream str_strm(str);
    while (str_strm >> tmp)
    {
        tokens.push_back(tmp);
    }
    return tokens;
}

// used to read our dataset/query files
void read_items(vector<Item> &data, const string &filename)
{
    string line;
    ifstream f;
    f.open(filename);
    if (f.is_open())
    {
        cout << "File " << filename << " successfully opened" << endl;
    }
    else
    {
        cout << "Error opening file " << filename << endl;
    }

    int count = 1;
    int catch_count = 0;
    int line_count = 1;
    while (!f.eof())
    {
        if (getline(f, line))
        {
            line_count++;
            vector<string> out;
            out = tokenize(line);
            Item item;
            item.id = out[0];
            item.null = false;

            // Μετατροπή των values του διανύσματος από string σε double
            for (size_t i = 1; i < out.size(); ++i)
            {
                try
                {
                    double xij_int = stod(out[i]);
                    item.xij.push_back(xij_int);
                }
                catch (...)
                {
                    cout << "[Line:" << line_count << "] Here is the catch: " << out[i] << "." << endl;
                    cout << "[Line:" << line_count << "] Type: " << typeid(out[i]).name() << endl;
                    catch_count++;
                }
            }
            data.push_back(item);
            count += 1;
        }

        /*
        else
        {
            cout << "Failed to read file at line " << count << endl;
        }
        */
    }
    // cout << "Catches: " << catch_count << endl;

    f.close();
    return;
}

// computes the Euclidean Distance between 2 Items of dimension d
double EuclideanDistance(const Item *p, const Item *q, const uint16_t &d)
{
    double sum = 0;
    for (int i = 0.0; i < d; i++)
    {
        sum += (p->xij[i] - q->xij[i]) * (p->xij[i] - q->xij[i]);
    }

    return sqrt(sum);
}

// used at sorting nearest neighbor points by distance
bool comparePairs(const std::pair<double, Item *> &x, const std::pair<double, Item *> &y)
{
    return (x.first < y.first);
}

// used to read data from config file
int extract_int_from_string(const string &str)
{
    stringstream ss;
    ss << str;
    string word;
    int temp_int;

    while (!ss.eof())
    {
        ss >> word;
        if (stringstream(word) >> temp_int)
        {
            ss.str(string());
            return temp_int;
        }
    }
    ss.str(string());
    return -1;
}

// used for testing purposes
template <typename K, typename V>
void print_map(std::unordered_map<K, V> const &m)
{
    for (auto const &pair : m)
    {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

// calculates Hamming distance between two bit-strings
int HammingDistance(unsigned int i1, unsigned int i2)
{
    int x = i1 ^ i2;
    int setBits = 0;

    while (x > 0)
    {
        setBits += x & 1;
        x >>= 1;
    }
    return setBits;
}

// used to calculate mean for the update(maximization) step of Lloyd's algorithm
// v1 works as an accumulator
vector<double> vector_mean(vector<double> &v1, const vector<double> &v2, int dimensions, int T)
{
    vector<double> result(dimensions, 0);
    for (int i = 0; i < dimensions; ++i)
    {
        result[i] = v1[i] + (v2[i] / T);
    }
    return result;
}

std::vector<std::pair<double, Item *>> brute_force_search(vector<Item> &dataset, Item *query, int N)
{
    int dimension = dataset[0].xij.size();

    // initialize a vector of N best candidates and distances represented as c++ pairs
    std::vector<std::pair<double, Item *>> knns;
    // Then initialize each pair with distance -> (max integer) and a null item
    for (int i = 0; i < N; i++)
        knns.push_back(std::make_pair(std::numeric_limits<double>::max(), new Item("null")));

    // For each item in dataset
    for (int j = 0; j < dataset.size(); j++)
    {
        // Calculate item's distance to the query item
        double distance = EuclideanDistance(query, &dataset[j], dimension);

        /*
        The last pair in the N-sized vector is the worst out of the N
        best candidates till now. If a better candidate is found,
        replace the last pair with the new one and re-sort the vector.
        */
        if (distance < knns[N - 1].first)
        {
            knns[N - 1].first = distance;
            if (knns[N - 1].second->null && knns[N - 1].second->id == "-1") // if it is a null item created just to initialize the N pairs of the vector.
                delete knns[N - 1].second;
            knns[N - 1].second = &dataset[j];
            std::sort(knns.begin(), knns.end(), comparePairs);
        }
    }
    return knns;
}

// Computes L2 norm between 2 points (x1,y1), (x2, y2)
double L2(double x1, double y1, double x2, double y2)
{
    double dist = pow((x1 - x2), 2.0) - pow((y1 - y2), 2.0);
    return sqrt(dist);
}

#endif
