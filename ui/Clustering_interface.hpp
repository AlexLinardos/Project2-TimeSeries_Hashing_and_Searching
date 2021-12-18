#ifndef CLUSTER_INTERFACE_HPP
#define CLUSTER_INTERFACE_HPP
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <algorithm>
#include "../includes/utils.hpp"

namespace Cli
{
    // class to hold parameters of Clustering
    class Cluster_params
    {
    private:
        // a map that will help us set the parameters from the command line
        std::map<std::string, std::string> param_set = {
            {"-i", "none"},
            {"-c", "none"},
            {"-o", "none"},
            {"-update", "none"},
            {"-assignment", "none"},
            {"-complete", "none"},
            {"-silhouette", "none"},
        };
        std::map<std::string, std::string>::iterator it;  // map iterator
        std::map<std::string, std::string>::iterator it2; // map iterator 2

        /* Makes sure the value of a given parameter in the configuration file can be converted to int and returns it.
            On a failure it returns default or -1 if there is no default. */
        int try_stoi(std::string param, int num)
        {
            try
            {
                int val = stoi(param);
                return val;
            }
            catch (...)
            {
                std::cout << "[Configuration ERROR] Given value " << param << " in configuration file cannot be converted to integer." << std::endl;
                if (num == 1)
                    return -1;
                if (num == 2)
                    return 3;
                if (num == 3)
                    return 4;
                if (num == 4)
                    return 10;
                if (num == 5)
                    return 3;
                if (num == 6)
                    return 2;

                return -2;
            }
        }
        

        // Sets parameter and confirms that all of the mandatory once have taken acceptable values.
        int set_and_confirm()
        {
            // iterate over map
            std::string curr_key;
            this->it = this->param_set.begin();
            while (it != this->param_set.end())
            {
                curr_key = it->first;
                if (it->second == "none")
                {
                    // handle optional
                    if ((curr_key != "-complete") && (curr_key != "-silhouette"))
                    {
                        std::cout << "[INPUT ERROR] Could not find value for mandatory parameter " << curr_key << ". If you entered this parameter then make sure you didn't forget anything before that." << std::endl;
                        return -1;
                    }
                }
                if (curr_key == "-i")
                    this->input_f = this->param_set.find("-i")->second;
                if (curr_key == "-c")
                    this->conf_f = this->param_set.find("-c")->second;
                if (curr_key == "-o")
                    this->output_f = this->param_set.find("-o")->second;
                if (curr_key == "-update")
                {
                    this->update = this->param_set.find("-update")->second;
                    std::string lc_update = lc(this->update);
                    if ((lc_update != "mean frechet") && (lc_update != "mean vector"))
                    {
                    	std::cout<<"Value "<<this->update<<" is not acceptable for parameter -update. Please enter Mean Frechet or Mean Vector."<<std::endl;
                    	return -1;
                    }
                }
                if (curr_key == "-assignment")
                {
                    this->assignment = this->param_set.find("-assignment")->second;
                    std::string lc_assignment = lc(this->assignment);
                    if ((lc_assignment != "classic") && (lc_assignment != "lsh") && (lc_assignment != "hypercube") && (lc_assignment != "lsh_frechet"))
                    {
                    	std::cout<<"Value "<<this->assignment<<" is not acceptable for parameter -assignment. Please enter Classic or LSH or Hypercube or LSH_Frechet."<<std::endl;
                    	return -1;
                    }
                }
                if (curr_key == "-complete")
                {
                    if (this->param_set.find("-complete")->second != "none")
                        this->complete = true;
                    else
                    	this->complete = false;
                }
                if (curr_key == "-silhouette")
                {
                    if (this->param_set.find("-silhouette")->second != "none")
                        this->silhouette = true;
                    else
                    	this->silhouette = false;
                }
                it++;
            }
            return 0;
        }

        // Reads parameter values from given configuration file
        int read_config()
        {
            std::string line;
            std::ifstream f;
            f.open(this->conf_f);
            if (f.is_open())
            {
                std::cout << "Configuration file " << this->conf_f << " successfully opened." << std::endl;
            }
            else
            {
                std::cout << "Error opening configuration file " << this->conf_f << std::endl;
                return -1;
            }

            int line_count = 0;
            while (!f.eof())
            {
                if (getline(f, line))
                {
                    line_count++;
                    std::vector<std::string> out;
                    out = tokenize(line);

                    if (line_count == 1)
                        this->clusters = this->try_stoi(out[1], 1);
                    if (line_count == 2)
                        this->L = this->try_stoi(out[1], 2);
                    if (line_count == 3)
                        this->k = this->try_stoi(out[1], 3);
                    if (line_count == 4)
                        this->M = this->try_stoi(out[1], 4);
                    if (line_count == 5)
                        this->d = this->try_stoi(out[1], 5);
                    if (line_count == 6)
                        this->probes = this->try_stoi(out[1], 6);
                }
            }
            // check if something went wrong with mandatory parameter K of K-medians
            if (this->clusters == -1)
            {
                std::cout << "[Configuration ERROR] Could not read mandatory parameter K of K-means." << endl;
                return -1;
            }

            f.close();
            return 0;
        }

    public:
        std::string input_f;    // input file
        std::string conf_f;     // configuration file
        int clusters;           // number of clusters aka K of K-medians
        int L;                  // number of hash tables for LSH
        int k;                  // number of vector hash functions for LSH
        int M;                  // max number of vectors to be checked for HC
        int d;                  // number of dimensions on which to project for HC
        int probes;             // max number of probes to be checked for HC
        std::string output_f;   // ouput file
        std::string update;     // method used for the update step
        std::string assignment; // method used for the assignment step
        bool complete;          // if given, prints clusters in detail in the output
        bool silhouette;        // if given, prints silhouette in the output
        bool success = true;    // to check if construction of object was 100% successful

        // Constructor
        Cluster_params(int argc, char *argv[])
        {
            if (argc == 1) // only for testing purposes
            {
                this->input_f = "datasets/nasd_input.csv";
                this->conf_f = "cluster.conf";
                this->output_f = "output.txt";
                this->update = "Mean Frechet"; // Mean Frechet - Mean Vector
                this->assignment = "Classic";  // Classic - LSH - Hypercube - LSH_Frechet
                this->complete = false;
                this->silhouette = false;
                if (this->read_config() == -1) // read the configuration file
                {
                    this->success = false;
                }
            }
            else
            {
                // pass user's parameter values into the map
                for (int i = 1; i < argc; i++)
                {
                    this->it = this->param_set.find(std::string(argv[i])); // try to find the parameter
                    if (this->it == this->param_set.end())                 // if NOT found in map
                    {
                        // it must be a value, so the previous argument should be the parameter or the word "Mean"
                        this->it2 = this->param_set.find(std::string(argv[i - 1]));
                        if (this->it2 == this->param_set.end()) // if not even the previous argument was an acceptable parameter
                        {
                            if ((std::string(argv[i-1])!="Mean") && (std::string(argv[i-1])!="mean")){
                            	std::cout << "[INPUT ERROR] Parameter " << std::string(argv[i]) << " is not compatible." << std::endl;
                            	this->success = false;
                            	break; // stop reading
                            }
                        }
                    }
                    else // if found pass the corresponding value to the map
                    {
                        // handle optional parameters
                        if ((this->it->first == "-complete") || (this->it->first == "-silhouette"))
                        {
                            this->it->second = "true";
                        }
                        else if (this->it->first == "-update")
                        {
                        	try 
                        	{
                        		this->it->second = std::string(argv[i+1])+" "+std::string(argv[i+2]);
                        	}
                        	catch(...)
                        	{
                        		std::cout<<"[INPUT ERROR] Failed to find value of -update parameter."<<std::endl;
                        		this->success = false;
                        	}
                        }
                        else
                        {
                        	try
                        	{
                        		this->it->second = std::string(argv[i + 1]);
                        	}
                        	catch(...)
                        	{
                        		std::cout<<"[INPUT ERROR] Failed to find value of -"<<std::string(argv[i])<<" parameter."<<std::endl;
                        		this->success = false;
                        	}
                            
                        }
                    }
                }
                if (this->set_and_confirm() == -1)
                {
                    this->success = false;
                }
                if (this->read_config() == -1) // read the configuration file
                {
                    this->success = false;
                }
            }
        }

        // Prints all parameters with their current values
        void print_NN_params()
        {
            std::cout << "-------------------------------------------" << std::endl
                      << "[Current parameters for Nearest Neighbor search]" << std::endl
                      << "Input file: " << this->input_f << std::endl
                      << "Configuration file: " << this->conf_f << std::endl
                      << "Output file: " << this->output_f << std::endl
                      << "Clusters: " << this->clusters << std::endl
                      << "L: " << this->L << std::endl
                      << "k: " << this->k << std::endl
                      << "M: " << this->M << std::endl
                      << "d': " << this->d << std::endl
                      << "probes: " << this->probes << std::endl
                      << "Update: " << this->update << std::endl
                      << "Assignment: " << this->assignment << std::endl
                      << "Complete: " << this->complete << std::endl
                      << "Silhouette: " << this->silhouette << std::endl
                      << "-------------------------------------------" << std::endl;
        }
    };
}
#endif
