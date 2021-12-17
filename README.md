# Project2 - Nearest Neighbour and Clustering for time series
### In this project we implement algorithms in C++ to perform Approximate Nearest Neighbour search and Clustering on datasets containing data series. For the implementation we experiment and use methods such as Locality Sensitive Hashing and Projection on Hypercube, as well as different distance metrics such as L2, Discrete Frechet and Continuous Frechet.

## Team members:
* Λινάρδος Αλέξανδρος - sdi1600093
* Αντώνης Μήλιος - sdi1600100

Github repository: https://github.com/AlexLinardos/Project2-TimeSeries_Hashing_and_Searching

## Repository organisation
### General
* __Makefile :__ Used to compile the programs with `make` command.
* __cluster.conf :__ Contains configuration parameters for clustering.
* __bin/ :__ Binary files are placed in this directory.
* __datasets/ :__ This directory contains the datasets used as input or as queries during the implementation.
* __outputs/ :__ Output files containing the results of the programs after a run are placed here.
 
### Code directories
For in-depth information about each file in these directories and the code please refer to the "In-depth Analysis" section further bellow.
* __TimeSeries-ANN/ :__ Contains code files necessary to implement Approximate Nearest Neighbour search with three different distance metrics (L2, Discrete Frechet, Continuous Frechet).
* __TimeSeries-Clustering/ :__
* __src/ :__ Contains the source (.cpp) files for both ANN search and Clustering.
* __interfaces/ :__ Contains header files for the implementation of the command line interfaces with which the user passes parameters and runs the programs.
* __misc/ :__ Contains header files with miscellaneous functions and classes used by our programs that don't have a direct connection to ANN or Clustering.

## Compilation and execution instructions
### Compilation options
* To compile the code for ANN search: `make search`
* To compile the code for Clustering: `make cluster`
To delete all the files in the bin directory use `make clean`

### Execution options
Please note that the parameters below noted with [M] are mandatory, while parameter "metric" noted with [M*] is mandatory only if LSH for curves (Frechet) algorithm option is selected. Non-mandatory paramateres will be set to default value if not given.
* To run ANN you are going to need the following parameters:
  1. input file[M] : the path to the input file (a dataset of curves), e.g. ../datasets/nasd_input.csv
  2. query file[M] : the path to the query file (a dataset of curves), e.g. ../datasets/nasd_query.csv
  3. k : the number of LSH hash functions to be used OR the dimensions of vectors d' for the projection on hypercube (dependent on algorithm)
  4. L : the number of LSH hash tables
  5. M : the maximum number of vectors to be checked while using Hypercube projection algorithm
  6. probes : the maximum number of hypercube vertices to be checked
  7. output file[M] : name of the file in which the results of the program will be written
  8. algorithm[M] : algorithm to be used for this run (LSH, Hypercube or LSH Frechet)
  9. metric[M*] : type of Frechet distance to be used if LSH for curves (LSH Frechet) algorithm is selected (discrete or continuous)
  10. delta : the size of the grid used in LSH for curves (LSH Frechet) algorithm
* To run Clustering you are going to need the following parameters:
  1. input file : the path to the input file (a dataset of curves), e.g. ./datasets/nasd_input.csv

To execute the __ANN search__ program follow this format:

`$./search –i <input file> –q <query file> –k <int> -L <int> -M <int> -probes
<int> -ο <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete
or continuous | only for –algorithm Frechet> -delta <double>`

To execute the __Clustering program__ follow this format:

`$./cluster –i <input file> –c <configuration file> -o <output file> -update <Mean
Frechet or Mean Vector> –assignment <Classic or LSH or Hypercube or LSH_Frechet>
-complete <optional> -silhouette <optional>`

__Examples :__

1. `./search -i ../datasets/nasd_input.csv -q ../datasets/nasd_query.csv -o output.txt -algorithm LSH`
2. `./search -i ../datasets/nasd_input.csv -q ../datasets/nasd_query.csv -o output.txt -algorithm Frechet -metric discrete -delta 2.0`
3. και ενα για cluster

## In-depth Analysis
In this section we will analyse our code file-by-file and talk about any possible implementation decisions taken during development.

#### TimeSeries-ANN/
#### TimeSeries-Clustering/
#### src/
#### interfaces/
#### misc/
