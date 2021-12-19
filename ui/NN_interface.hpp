#ifndef NN_INTERFACE_HPP
#define NN_INTERFACE_HPP
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include "../includes/utils.hpp"

namespace NNi
{
    // class to hold parameters of Nearest Neighbor search
    // values are initialized to default
    class NN_params
    {
    private:
        // a map that will help us set the parameters from the command line
        std::map<std::string, std::string> param_set = {
            {"-i", "none"},
            {"-q", "none"},
            {"-k", "none"},
            {"-L", "none"},
            {"-M", "none"},
            {"-probes", "none"},
            {"-o", "none"},
            {"-algorithm", "none"},
            {"-metric", "none"},
            {"-delta", "none"}};
        std::map<std::string, std::string>::iterator it; // map iterator

        /* Makes sure the value of a given parameter can be converted to int and sets it.
            On a failure it marks the value as "none" so as to be set to default later if able. */
        void try_stoi(std::string param)
        {
            try
            {
                int val = stoi(this->param_set.find(param)->second);
                if (param == "-k")
                    this->k = val;
                else if (param == "-L")
                    this->L = val;
                else if (param == "-M")
                    this->M = val;
                else if (param == "-probes")
                    this->probes = val;
            }
            catch (...)
            {
                std::cout << "[INPUT ERROR] Given value " << this->param_set.find(param)->second << " of parameter " << param << " cannot be converted to integer." << std::endl;
                this->param_set.find(param)->second = "none";
            }
        }

        /* Makes sure the value of a given parameter can be converted to double and sets.
            On a failure it resets the value to "none". */
        void try_stod(std::string param)
        {
            try
            {
                double val = stod(this->param_set.find(param)->second);
                if (param == "-delta")
                    this->delta = val;
            }
            catch (...)
            {
                std::cout << "[INPUT ERROR] Given value " << this->param_set.find(param)->second << " of parameter " << param << " cannot be converted to double." << std::endl;
                this->param_set.find(param)->second = "none";
            }
        }

        /* Sets default value to corresponding command line parameter.
            Default value of parameter k changes according to algorithm.
            Prints error message on failure. */
        void set_default(std::string param, std::string algo)
        {
            if (param == "-k")
            {
                if (algo == "LSH")
                {
                    this->k = 4;
                    this->param_set.find(param)->second = "4";
                }
                else
                {
                    this->k = 14;
                    this->param_set.find(param)->second = "14";
                }
            }
            else if (param == "-L")
            {
                this->L = 5;
                this->param_set.find(param)->second = "5";
            }
            else if (param == "-M")
            {
                this->M = 10;
                this->param_set.find(param)->second = "10";
            }
            else if (param == "-probes")
            {
                this->probes = 2;
                this->param_set.find(param)->second = "2";
            }
            else
            {
                std::cout << "Parameter " << param << " has no default value." << std::endl;
            }
        }

        /* Sets parameter and confirms that all of them have taken values.
            If not, it places default values when able.
            Returns -1 on failure to set needed parameters.*/
        int set_and_confirm()
        {
            // first of all set the algorithm and the metric
            std::string algo = this->param_set.find("-algorithm")->second;
            std::string algo_lc = lc(algo);
            if ((algo_lc != "lsh") && (algo_lc != "hypercube") && (algo_lc != "frechet"))
            {
                std::cout << "[INPUT ERROR] Algorithm " << algo << " is not recognisable." << std::endl;
                return -1;
            }
            else
            {
                this->algorithm = algo; // set algorithm value
                if (algo_lc == "frechet")
                {
                    std::string met = this->param_set.find("-metric")->second;
                    std::string met_lc = lc(met);
                    if ((met_lc != "discrete") && (met_lc != "continuous"))
                    {
                        std::cout << "[INPUT ERROR] Metric " << met << " is not recognisable." << std::endl;
                        return -1;
                    }
                    else
                    {
                        this->metric = met;
                    }
                }
            }

            // iterate over map
            std::string curr_key;
            this->it = this->param_set.begin();
            while (it != this->param_set.end())
            {
                curr_key = it->first;
                // if a parameter has no value see if it can be set to a default
                if (it->second == "none")
                {
                    if ((curr_key == "-k") || (curr_key == "-L") || (curr_key == "-M") || (curr_key == "-probes"))
                    {
                        this->set_default(curr_key, algo);
                        std::cout << "Parameter " << curr_key << " has been set to default value." << std::endl;
                    }
                    else
                    {
                        std::cout << "[INPUT ERROR] Could not find value for mandatory parameter " << curr_key << "." << std::endl;
                        return -1;
                    }
                }
                // else set the value while doing type conversion if needed
                else
                {
                    if (curr_key == "-k")
                        this->try_stoi("-k");
                    else if (curr_key == "-L")
                        this->try_stoi("-L");
                    else if (curr_key == "-M")
                        this->try_stoi("-M");
                    else if (curr_key == "-probes")
                        this->try_stoi("-probes");
                    else if (curr_key == "-delta")
                        this->try_stod("-delta");
                    else if (curr_key == "-i")
                        this->input_f = this->param_set.find("-i")->second;
                    else if (curr_key == "-q")
                        this->query_f = this->param_set.find("-q")->second;
                    else if (curr_key == "-o")
                        this->output_f = this->param_set.find("-o")->second;
                }
                it++;
            }

            // final check for potential stoi fails (we should set to default)
            this->it = this->param_set.begin();
            while (it != this->param_set.end())
            {
                curr_key = it->first;
                // if a parameter has no value see if it can be set to a default
                if (it->second == "none")
                {
                    if ((curr_key == "-k") || (curr_key == "-L") || (curr_key == "-M") || (curr_key == "-probes"))
                    {
                        this->set_default(curr_key, algo);
                        std::cout << "Parameter " << curr_key << " has been set to default value." << std::endl;
                    }
                    else if (curr_key == "-delta")
                    {
                        this->param_set.find("-delta")->second = "0.0";
                    }
                    else
                    {
                        std::cout << "[INPUT ERROR] Could not find value for mandatory parameter " << curr_key << "." << std::endl;
                        return -1;
                    }
                }
                it++;
            }
            return 0;
        }

    public:
        std::string input_f;   // input file
        std::string query_f;   // query file
        int k;                 // number of hi hash functions for LSH or d' dimensions for HC
        int L;                 // number of hash tables for LSH
        int M;                 // max number of vectors to be checked for HC
        int probes;            // max number of probes to be checked for HC
        std::string output_f;  // ouput file
        std::string algorithm; // algorithm to be used for this run
        std::string metric;    // metric to be used for this run (only for -algorithm Frechet)
        double delta;          // delta parameter for LSH of curves
        bool success = true;   // to check if construction of object was 100% successful

        // Constructor
        NN_params(int argc, char *argv[])
        {
            if (argc == 1) // only for testing purposes
            {
                this->input_f = "datasets/nasd_input.csv";
                this->query_f = "datasets/nasd_query.csv";
                this->k = 4;
                this->L = 5;
                this->M = 0;
                this->probes = 0;
                this->output_f = "outputs/output.txt";
                this->algorithm = "Frechet"; // LSH - Hypercube - Frechet
<<<<<<< HEAD
                this->metric = "continuous";       // discrete - continuous
=======
                this->metric = "discrete";   // discrete - continuous
>>>>>>> acc225e10fafd8c683392d131c086c201c8f4458
                this->delta = 0.0;
            }
            else if (argc % 2 == 1) // normal run
            {
                // pass user's parameter values into the map
                for (int i = 1; i < argc; i = i + 2)
                {
                    this->it = this->param_set.find(std::string(argv[i])); // try to find the parameter
                    if (this->it == this->param_set.end())                 // if user entered a parameter that does not exist
                    {
                        std::cout << "[INPUT ERROR] Parameter " << argv[i] << " is not compatible." << std::endl;
                        this->success = false;
                        break; // stop reading
                    }
                    else // if found pass the corresponding value to the map
                    {
                        this->it->second = std::string(argv[i + 1]);
                    }
                }

                // set parameters to corresponding values of map (or default if not mandatory)
                int result = this->set_and_confirm();
                if (result == -1)
                    this->success = false;
            }
            else
            {
                std::cout << "[INPUT ERROR] Missing command line parameter." << std::endl;
                this->success = false;
            }
        }

        // Prints all parameters with their current values
        void print_NN_params()
        {
            std::cout << "-------------------------------------------" << std::endl
                      << "[Current parameters for Nearest Neighbor search]" << std::endl
                      << "Input file: " << this->input_f << std::endl
                      << "Query file: " << this->query_f << std::endl
                      << "Output file: " << this->output_f << std::endl
                      << "Algorithm: " << this->algorithm << std::endl
                      << "Metric: " << this->metric << std::endl
                      << "k: " << this->k << std::endl
                      << "L: " << this->L << std::endl
                      << "M: " << this->M << std::endl
                      << "probes: " << this->probes << std::endl
                      << "delta: " << this->delta << std::endl
                      << "-------------------------------------------" << std::endl;
        }
    };

}

#endif
