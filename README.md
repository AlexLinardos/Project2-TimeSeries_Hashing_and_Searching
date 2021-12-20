# Project2 - Nearest Neighbour and Clustering for time series
### In this project we implement algorithms in C++ to perform Approximate Nearest Neighbour search and Clustering on datasets containing time series. For the implementation we experiment and use methods such as Locality Sensitive Hashing and Projection on Hypercube, as well as different distance metrics such as L2, Discrete Frechet and Continuous Frechet.

## Team members:
* Λινάρδος Αλέξανδρος - sdi1600093
* Αντώνης Μήλιος - sdi1600100

Github repository: https://github.com/AlexLinardos/Project2-TimeSeries_Hashing_and_Searching

## Table of contents
* [Repository organisation](#repository-organisation)
* [Compilation and execution instructions](#compilation-and-execution-instructions)
* [In-depth Analysis](#in-depth-analysis)
* [Optimization experiments and parameter tuning](#optimization-experiments-and-parameter-tuning)
* [Evaluation and final thoughts](#evaluation-and-final-thoughts)

## Repository organisation
### General
* __Makefile :__ Used to compile the programs with `make` command.
* __cluster.conf :__ Contains configuration parameters for clustering. When changing the values make sure you *__leave a space before and after the given value__*.
* __bin/ :__ Binary files are placed in this directory.
* __datasets/ :__ This directory contains the datasets used as input or as queries during the implementation.
* __outputs/ :__ This directory contains output files with the results of some program runs.
 
### Code directories
*For in-depth information about each file in these directories and the code please refer to the [In-depth Analysis](#in-depth-analysis) section further bellow.*
* __TimeSeries-ANN/ :__ Contains code files necessary to implement Approximate Nearest Neighbour search with three different distance metrics (L2, Discrete Frechet, Continuous Frechet).
* __TimeSeries-Clustering/ :__
* __src/ :__ Contains the source (.cpp) files for both ANN search and Clustering.
* __ui/ :__ Contains header files for the implementation of the command line interfaces with which the user passes parameters and runs the programs.
* __includes/ :__ Contains header files with various utility functions and classes used by our programs that don't have a direct connection to ANN or Clustering.

## Compilation and execution instructions
### Compilation options
* To compile the code for ANN search: `make search`
* To compile the code for Clustering: `make cluster`
* To compile everything: `make` 
<br>To delete all the files in both the "bin" and the "outputs" directories use `make clean`
<br>To delete all the files in the "bin" directory use `make clean_bin`

### Execution options
Please note that the parameters below noted with [__M__] are mandatory, while parameter "metric" noted with [__M*__] is mandatory only if LSH for curves (Frechet) algorithm option is selected. Non-mandatory paramateres will be set to default value if not given.
* To run ANN you are going to need the following parameters:
  1. input file[__M__] : the path to the input file (a dataset of curves), e.g. ../datasets/nasd_input.csv
  2. query file[__M__] : the path to the query file (a dataset of curves), e.g. ../datasets/nasd_query.csv
  3. k : the number of LSH hash functions to be used OR the dimensions of vectors (aka d') for the projection on hypercube (dependent on algorithm)
  4. L : the number of LSH hash tables
  5. M : the maximum number of vectors to be checked while using Hypercube projection algorithm
  6. probes : the maximum number of hypercube vertices to be checked
  7. output file[__M__] : the path to the file in which the results of the program will be written
  8. algorithm[__M__] : algorithm to be used for this run (LSH, Hypercube or LSH Frechet)
  9. metric[__M*__] : type of Frechet distance to be used if LSH for curves (LSH Frechet) algorithm is selected (discrete or continuous)
  10. delta : the size of the grid used in LSH for curves (LSH Frechet) algorithm
* To run Clustering you are going to need the following parameters:
  1. input file[__M__] : the path to the input file (a dataset of curves), e.g. ./datasets/nasd_input.csv
  2. configuration file[__M__] : the path to the configuration (.conf) file, e.g. ./cluster.conf <br>
      A configuration file should include the following parameters and in this exact order:
      * number_of_clusters[__M__] : number of clusters to be created (aka K of K-medians)
      * number_of_vector_hash_tables : the number of LSH hash tables
      * number_of_vector_hash_functions : the number of LSH hash functions to be used
      * max_number_M_of_hypercube : the maximum number of vectors to be checked
      * number_of_hypercube_dimensions : the dimensions of vectors (aka d') for the projection on hypercube
      * number_of_probes : the maximum number of hypercube vertices to be checked
  3. output file[__M__] : the path to the file in which the results of the program will be written
  4. update[__M__] : method to be used during the update step of clustering (Mean Frechet or Mean Vector)
  5. assignment[M] : method to be used during the assignment step of clustering (Classic, LSH, Hypercube or LSH_Frechet)
  6. complete: if given, prints clusters in detail in the output
  7. silhouette: if given, prints silhouette in the output

To execute the __ANN search__ program follow this format:

`$./search –i <input file> –q <query file> –k <int> -L <int> -M <int> -probes
<int> -ο <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete
or continuous | only for –algorithm Frechet> -delta <double>`

To execute the __Clustering program__ follow this format:

`$./cluster –i <input file> –c <configuration file> -o <output file> -update <Mean
Frechet or Mean Vector> –assignment <Classic or LSH or Hypercube or LSH_Frechet>
-complete <optional> -silhouette <optional>`

__Examples :__

1. ./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -k 1 -L 1 -o outputs/output.txt -algorithm LSH
2. `./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -o outputs/output.txt -algorithm Frechet -metric discrete -delta 2.0`
3. `./bin/cluster -i datasets/nasd_input.csv -c cluster.conf -o outputs/output.txt -update Mean Vector -assignment Classic`
4. `./bin/cluster -i datasets/nasd_input.csv -c cluster.conf -o outputs/output.txt -update Mean Frechet -assignment LSH_Frechet -complete -silhouette`

## In-depth Analysis
In this section we will analyse our code file-by-file and talk about any possible implementation decisions taken during development.

### TimeSeries-ANN/
* #### TimeSeries-ANN/L2
  This directory includes header files that implement the LSH and Hypercube projection methods for approximate nearest neighbour searching using euclidean distance as metric. They were implemented as part of a previous project so for more info please refer to the following GitHub repository: https://github.com/AlexLinardos/Project1-LSH-and-Clustering
* #### TimeSeries-ANN/DiscreteFrechet
  1. __disc_Frechet.hpp__ : This header file contains a function -called discrete_frechet()- that *__computes the discrete Frechet distance__* between two curves that are given as parameters. The computation is done with dynamic programming via the function compute_c() that is included in the same file. It is important to note that discrete_frechet() returns the whole table of values that is computed with the dynamic programming approach and not just the final result. Also included in this file is function search_exactNN() that -given a query curve and a dataset- *__searches for the exact nearest neighbour of the query inside the dataset__*. It returns a pair that contains a pointer to the nearest neighbour curve found and the Frechet distance from the query. Note that this is a different function from the one with the same name that is included in the ContinuousFrechet directory as this one uses *discrete* Frechet distance as a metric (while the other uses continuous Frechet distance).
  2. __discF_LSH.hpp__ : This header file contains the LSH class (see comments on code for explanation of each parameter of constructor). Upon creation, an object of this class has everything needed to *__hash and append each curve of a given dataset into multiple hash tables__* via the dataset_hashing() method. Note that the hash tables are filled with Association objects, which are triplets of a *curve*, its *grid curve* (produced by the produce_h() method) and the grid curve as a 1d *vector* (aka a flattened version that is produced by the concat_points() method). The produce_h() method uses random shifts distributed uniformly in [0, δ) to produce different grids (as many as the hash tables) and *__snap the curves__* onto them using the following formula: xi' = floor((x-t)/δ + 1/2)δ + t, yi' = floor((x-t)/δ + 1/2)δ + t, where δ is delta parameter and t is a randomized shift each time. After this procedure, the associations are stored using their vector as key into the 1d tables with the help of LSH for vectors. Finally, the search_ANN() method hashes a given query curve in the same way and then *__finds an approximate nearest neighbour__* of the query inside the dataset by searching the curves that are in the same hash bucket as itself and in all hash tables. The comparisons are done using discrete Frechet distance and the method returns a pair of a pointer to the approximately nearest curve and the discrete Frechet distance from it.<br><br>
  __Notable implementation decisions__:
  * Function for discrete Frechet distance calculation returns the whole table of the dynamic programming approach instead of just the final value.
  * Shifts t are different for each dimension.
  * After snapping curves to a grid, we replace lost points by applying padding so as all vectors have the same length before using LSH to store them in hash tables.<br><br>
  
* #### TimeSeries-ANN/ContinuousFrechet
  1. __cont_Frechet.hpp__ : This header file contains functions that were used to help us use the code in the "ContinuousFrechet/Fred" directory while treating it as a __"black box"__. They (a) convert our objects to objects of the classes that are implemented in the "ContinuousFrechet/Fred" directory and (b) use functions that are implemented in the "ContinuousFrechet/Fred" directory to calculate the Continuous Frechet Distance. Also included in this file is the function search_exactNN() that -given a query curve and a dataset- *__searches for the exact nearest neighbour of the query inside the dataset by using the continuous Frechet distance__*.
  2. __contF_LSH.hpp__ : This header file contains the LSH class (see comments on code for explanation of each parameter of constructor). Upon creation, an object of this class has everything needed to *__hash and append each curve of a given dataset into multiple hash tables__* via the dataset_hashing() method. Note that the hash tables are filled with Association objects, which are triplets of a *curve*, its *grid curve* (produced by the snap_to_1dgrid() method) and the grid curve as a 1d *vector* (aka a flattened version). The snap_to_1dgrid() method uses random shifts distributed uniformly in [0, δ) to produce different grids (as many as the hash tables) and *__snap the curves__* onto them using the following formula: xi' = floor((x-t)/δ + 1/2)δ + t, yi' = floor((x-t)/δ + 1/2)δ + t, where δ is delta parameter and t is a randomized shift each time. After this procedure, the associations are stored using their vector as key into the 1d tables with the help of LSH for vectors. Finally, the search_ANN() method hashes a given query curve in the same way and then *__finds an approximate nearest neighbour__* of the query inside the dataset by searching the curves that are in the same hash bucket as itself and in all hash tables. The comparisons are done using continuous Frechet distance and the method returns a pair of a pointer to the approximately nearest curve and the continuous Frechet distance from it.
  3. __Fred/__ : The code contained in this directory was given to us to help us with the calculation of the continuous Frechet distance and was used as a "black box". *Source*: https://github.com/derohde/Fred/

### TimeSeries-Clustering/
  1. __curve_clustering.hpp__ : This header file contains a Clustering class. Objects of this class contain all the required methods/algorithms for centroid initialization, assignment and update, thus enabling us to perform clustering on a dataset of curves. Note that all of the methods in this class use the discrete Frechet distance metric, since they are made for curve clustering. The initialize_pp() method performs the *__k-means++ initialization process to initialize the centroids of a given number of clusters__*. For that, it uses the methods in the private section to help with the calculations of the distances from nearest centroids and the probability of each point to be chosen as next centroid. For *__the assignment step__* (where we assign each item to a cluster) there are two options : the classic assignment (Lloyds_assignment()) method and the assignment by reverse LSH range search (Range_dfLSH_assignment()) method. *__The update step__* of the clustering is done via the update_centers() method. Finally, *__combining all of the steps of clustering to implement a complete clustering method__* are the Classic_Clustering() and the Reverse_Assignment_LSH_Clustering() methods. These methods differ on the assignment method that they use, where the first uses classic assignment and the second uses the assignment by reverse LSH range search. Both of these methods have a maximum iteration threshold taken as a parameter to avoid a case of an infinite loop if the centroids don't converge. 
  2. __vector_clustering.hpp__ : This header file contains a Clustering class. Objects of this class contain all the required methods/algorithms for centroid initialization, assignment and update, thus enabling us to perform clustering on a dataset of vectors. Note that all of the methods in this class use the L2 distance metric, since they are made for vector clustering. The initialize_pp() method performs the *__k-means++ initialization process to initialize the centroids of a given number of clusters__*. For that, it uses the methods in the private section to help with the calculations of the distances from nearest centroids and the probability of each point to be chosen as next centroid. For *__the assignment step__* (where we assign each item to a cluster) there are three options : the classic assignment (Lloyds_assignment()) method, the assignment by reverse LSH range search (Range_LSH_assignment()) method and the assignment by reverse Hypercube range search method. *__The update step__* of the clustering is done via the update_centers() method. Finally, *__combining all of the steps of clustering to implement a complete clustering method__* are the Classic_Clustering(), Reverse_Assignment_LSH_Clustering() and Reverse_Assignment_HC_Clustering() methods. These methods differ on the assignment method that they use, where the first uses classic assignment, the second uses the assignment by reverse LSH range search and the third uses the assignment by reverse Hypercube range search. Both of these methods have a maximum iteration threshold taken as a parameter to avoid a case of an infinite loop if the centroids don't converge.
  3. __silhouette.hpp__ : 


### src/
  1. __search_main.cpp__ : Main code file that executes the nearest neighbour searching process according to the command line parameters given by the user. It checks the value of the -algorithm parameter and executes the respective algorithm. If the algorithm is set to "Frechet" it also check the -metric parameter and executes accordingly. At the end, it produces the output file with the results of the run and also prints evaluation metrics on std out.

### ui/
  1. __NN_interface.hpp__: This header file contains the class NN_params. After construction, an object of this class has all the necessary attributes and methods to implement a fully-fledged command line interface for the ANN search program that includes methods for reading a command, checking if it is valid, finding parameters, confirming that their values are acceptable, storing them and passing default values into possible non-mandatory parameters that were not given. After its construction, all parameters needed for the ANN search algorithm will be stored in its attributes.
  2. __Clustering_interface.hpp__: This header file contains the class Cluster_params. After construction, an object of this class has all the necessary attributes and methods to implement a fully-fledged command line interface for the Clustering program that includes methods for reading a command, checking if it is valid, finding parameters, confirming that their values are acceptable, storing them and passing default values into possible non-mandatory parameters that were not given. After its construction, all parameters needed for the ANN search algorithm will be stored in its attributes.

### includes/
  1. __curves.hpp__: This header file contains a class that is used throughout our code to make it more intuitive. That is the Curve2d class. It helps us create objects that *__simulate curves in the 2-dimensional space__* by having an id, a vector of 2-dimensional points and some attributes that help us perform clustering (such as the "cluster" attribute that is the index of the cluster to which said curve is assigned). The 2-dimensional points are constructed by the class just above, named Point2d. This is again a very simple class that is there only to make the rest of our code more intuitive. The point2d_L2() function calculates the euclidean distance between two of those points. Finally, included here is the delta_tuning() function that is responsible for *__calculating the delta parameter__* in case it is not given by the user and the identical_curves() function to *__help us implement a quering trick for LSH__*. For more info about both of these function please refer to the [Optimization experiments and parameter tuning](#optimization-experiments-and-parameter-tuning) section.
  2. __mean_curve.hpp__: This header file contains a function mean_curve() *__to compute the mean curve of two given curves__*. The mean curve is filtered in order to be simplified. The computation of a mean curve requires an *__optimal traversal__* of the given curves, which is computed by the function optimal_traversal() in the same file. Finally, the mean_of_curves() function returns *__the mean curve  of multiple given curves__*. <br><br>
  __Notable implementation decisions__:
  * The algorithm to compute the optimal traversal returns it in reverse (because vectors do not and should not have a way to push an element to the front). We decided not to reverse the returned traversal in order to not add additional computational weight to the program.
  * The mean_of_curves() function computes the needed mean curve by simulating a binary tree. At first, it begins with the last two curves given in the std::vector that is given as parameter and calculates their mean. Then that mean curve is stored in an std::vector and we decrease the index showing our remaining given curves by 2. This procedure repeats until either (a) the "remaining" index reaches 0 -meaning we have calculated half the mean curves of each concecutive pair- or (b) the "remaining" index reaches 1 -meaning the number of curves is odd, so one of them does not have a pair-. In the second case, we push the remaining curve together with the previously computed mean curves. After this initial "loading" of mean curves in the aforementioned std::vector, we repeat the following process until there is only 1 remaining curve:
    * Push the mean curve of each concecutive pair of curves in the std::vector
    * Decrease the "remaining" index by 2 so as to not select newly added mean curves at the same iteration during which they were pushed in the std::vector
    * The pair of curves of which we calculated the mean is removed.<br>
  At the end of this procedure, the only remaining curve will be the mean curve of all the given curves.

## Optimization experiments and parameter tuning

__Quering trick__ (for Frechet methods): When we find the hash bucket of the query, we check the bucket for the existance of identical grid-curves so as to skip the calculation of the Frechet distance for the non-identical ones. We implemented this trick to speed up the procedure but there was no noticable improvement (at least for our dataset) for the Discrete Frechet ANN approach as the calculation of the metric is fast enough. That said, the current implementation of the Discrete Frechet ANN runs without the querying trick but there is always the option to use it by calling the search_ANN() method with the "quering_trick" boolean parameter set to "true" (from search_main.cpp). 
We believe that this trick can be of help, at the Continuous Frechet ANN approach in order to avoid the demanding calculation of the metric, but testing and evaluating this enough to have concrete proof was hard due to the continuous frechet calculation being slow. That said, the current implementation of the Continuous Frechet ANN runs with the querying trick on.

__Padding__ : To be able to use the LSH for vectors method we had to apply padding (meaning to replace any points that were cut) to the grid-curves in order to have all the vectors be of the same length. We chose to use the value 10000 for the padding. If you wish to use another dataset keep in mind that it would be best to replace that value with a value big enough to be outside of the range of values in your dataset.

__Delta tuning__: When not given by the user, this parameter is computed automatically by the program according to the given data so as to approximate the mean distance between two concective vertices in the curves of the dataset. In order to achieve that, the delta_tuning() function (located in curves.hpp), for a number of iterations, selects random curve pairs, sums the mean distances of their vertices and calculates the average value for all pairs.

### Best parameters for each algorithm
After experimenting we concluded that the best parameters for each algorithm are as follows:<br>
* *Vector LSH* : L=1, k=1, tablesize=datasetsize/2,  w=average L2 distance between dataset curves (as vectors)
* *Hypercube*: k=1. M=30, probes=1, w=average L2 distance between dataset curves (as vectors)
* *discrete Frechet LSH* : L=6, delta=avg dist between curve vertices, tablesize=dataset_size/8, threshold=dataset_size/4, querying trick=false
* *continuous Frechet LSH* : L=1, delta=avg dist between curve vertices, tablesize=dataset_size/4, threshold=dataset_size/4, querying trick=true

## Evaluation and final thoughts

### Nearest Neighbour search
For the nearest neighbour search we used the following metrics to evaluate the performance of our code:
* Maximum Approximation Factor (MAF)
* average(distanceApproximate / distanceTrue) (printed in std out)
* timeApproximateAverage/timeTrueAverage (printed in std out)

Performance tests per method:<br>
* __LSH/Hypercube for vectors__ : As expected, these methods did not perform as good as they did in our previous project. While they are not unusable (they finish slightly faster than the brute force method without being too far off), it is pretty clear that the complexity of the curves requires a different approach in order to have satisfying results.
* __LSH for curves with discrete Frechet__ : Definitely an upgrade compare to the previous method and certainly the best overall. This method performs a lot faster than the brute force method and reliably returns great approximations.
* __LSH for curves with continuous Frechet__ : While incredibly accurate, this method is extremely slow. Its low speed maybe caused by the complexity of the metric, as well as the fact that we were indicated to use "black box" functions to call the already implemented code that was provided to us (Fred folder). Even after a lot of filtering to reduce the complexity of the curves, the algorithm remains pretty slow.

./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -k 1 -L 1 -o outputs/output.txt -algorithm LSH




