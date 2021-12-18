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
* __outputs/ :__ This directory contains output files with the results of some program runs.
 
### Code directories
*For in-depth information about each file in these directories and the code please refer to the "In-depth Analysis" section further bellow.*
* __TimeSeries-ANN/ :__ Contains code files necessary to implement Approximate Nearest Neighbour search with three different distance metrics (L2, Discrete Frechet, Continuous Frechet).
* __TimeSeries-Clustering/ :__
* __src/ :__ Contains the source (.cpp) files for both ANN search and Clustering.
* __ui/ :__ Contains header files for the implementation of the command line interfaces with which the user passes parameters and runs the programs.
* __includes/ :__ Contains header files with various utility functions and classes used by our programs that don't have a direct connection to ANN or Clustering.

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
  7. output file[M] : the path to the file in which the results of the program will be written
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

1. `./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -o outputs/output.txt -algorithm LSH`
2. `./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -o outputs/output.txt -algorithm Frechet -metric discrete -delta 2.0`
3. και ενα για cluster

## In-depth Analysis
In this section we will analyse our code file-by-file and talk about any possible implementation decisions taken during development.

### TimeSeries-ANN/
* #### TimeSeries-ANN/L2
  This directory includes header files that implement the LSH and Hypercube projection methods for approximate nearest neighbour searching using euclidean distance as metric. They were implemented as part of a previous project so for more info please refer to the following GitHub repository: https://github.com/AlexLinardos/Project1-LSH-and-Clustering
* #### TimeSeries-ANN/DiscreteFrechet
  1. __disc_Frechet.hpp__ : This header file contains a function -called discrete_frechet()- that *__computes the discrete Frechet distance__* between two curves that are given as parameters. The computation is done with dynamic programming via the function compute_c() that is included in the same file. It is important to note that discrete_frechet() returns the whole table of values that is computed with the dynamic programming approach and not just the final result. Also included in this file is function search_exactNN() that -given a query curve and a dataset- *__searches for the exact nearest neighbour of the query inside the dataset__*. It returns a pair that contains a pointer to the nearest neighbour curve found and the Frechet distance from the query. Note that this is a different function from the one with the same name that is included in the ContinuousFrechet directory as this one uses *discrete* Frechet distance as a metric (while the other uses continuous Frechet distance).
  2. __discF_LSH.hpp__ : This header file contains the LSH class (see comments on code for explanation of each parameter of constructor). Upon creation, an object of this class has everything needed to *__hash and append each curve of a given dataset into multiple hash tables__* via the dataset_hashing() method. Note that the hash tables are filled with Association objects, which are triplets of a *curve*, its *grid curve* (produced by the produce_h() method) and the grid curve as a 1d *vector* (aka a flattened version that is produced by the concat_points() method). The produce_h() method uses random shifts distributed uniformly in [0, δ) to produce different grids (as many as the hash tables) and *__snap the curves__* onto them using the following formula: xi' = floor((x-t)/δ + 1/2)δ + t, yi' = floor((x-t)/δ + 1/2)δ + t, where δ is delta parameter and t is a randomized shift each time. After this procedure, the associations are stored by their vector into the 1d tables with the help of LSH for vectors. Finally, the search_ANN() method hashes a given query curve in the same way and then *__finds an approximate nearest neighbour__* of the query inside the dataset by searching the curves that are in the same hash bucket as itself and in all hash tables. The comparisons are done using discrete Frechet distance and the method returns a pair of a pointer to the approximately nearest curve and the discrete Frechet distance from it.<br><br>
  __Notable implementation decisions__:
  * Function for discrete Frechet distance calculation returns the whole table of the dynamic programming approach instead of just the final value.
  * Shifts t are different for each dimension.
  * After snapping curves to a grid, we replace lost points by applying padding so as all vectors have the same length before using LSH to store them in hash tables.<br><br>
  
* #### TimeSeries-ANN/ContinuousFrechet
  1. ...

### TimeSeries-Clustering/

### src/

### ui/
  1. __NN_interface.hpp__: This header file contains the class NN_params. After construction, an object of this class has all the necessary attributes and methods to implement a fully-fledged command line interface for the ANN search program that includes methods for reading a command, checking if it is valid, finding parameters, confirming that their values are acceptable, storing them and passing default values into possible non-mandatory parameters that were not given. After its construction, all parameters needed for the ANN search algorithm will be stored in its attributes.
  2. __Clustering_interface.hpp__: ... 

### includes/
  1. __curves.hpp__: This header file contains a class that is used throughout our code to make it more intuitive. That is the Curve2d class. It helps us create objects that *__simulate curves in the 2-dimensional space__* by having an id and a vector of 2-dimensional points. The 2-dimensional points are constructed by the class just above, named Point2d. This is again a very simple class that is there only to make the rest of our code more intuitive. Finally, the point2d_L2() function calculates the euclidean distance between two of those points.
  2. __mean_curve.hpp__: This header file contains a function mean_curve() *__to compute the mean curve of two given curves__*. The mean curve is filtered in order to be simplified. The computation of a mean curve requires an *__optimal traversal__* of the given curves, which is computed by the function optimal_traversal() in the same file. Finally, the mean_of_curves() function returns *__the mean curve  of multiple given curves__*. <br><br>
  __Notable implementation decisions__:
  * The algorithm to compute the optimal traversal returns it in reverse (because vectors do not and should not have a way to push an element to the front). We decided not to reverse the returned traversal in order to not add additional computational weight to the program.
  * εδώ για το mean_of_curves
