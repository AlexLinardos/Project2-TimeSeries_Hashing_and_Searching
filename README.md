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
* [Final thoughts](#final-thoughts)

## Repository organisation
### General
* __Makefile :__ Used to compile the programs with `make` command.
* __cluster.conf :__ Contains configuration parameters for clustering. When changing the values make sure you *__leave a space before and after the given value__*.
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

1. `./bin/search -i datasets/nasd_input.csv -q datasets/nasd_query.csv -o outputs/output.txt -algorithm LSH`
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
  2. __Clustering_interface.hpp__: This header file contains the class Cluster_params. After construction, an object of this class has all the necessary attributes and methods to implement a fully-fledged command line interface for the Clustering program that includes methods for reading a command, checking if it is valid, finding parameters, confirming that their values are acceptable, storing them and passing default values into possible non-mandatory parameters that were not given. After its construction, all parameters needed for the ANN search algorithm will be stored in its attributes.

### includes/
  1. __curves.hpp__: This header file contains a class that is used throughout our code to make it more intuitive. That is the Curve2d class. It helps us create objects that *__simulate curves in the 2-dimensional space__* by having an id and a vector of 2-dimensional points. The 2-dimensional points are constructed by the class just above, named Point2d. This is again a very simple class that is there only to make the rest of our code more intuitive. Finally, the point2d_L2() function calculates the euclidean distance between two of those points.
  2. __mean_curve.hpp__: This header file contains a function mean_curve() *__to compute the mean curve of two given curves__*. The mean curve is filtered in order to be simplified. The computation of a mean curve requires an *__optimal traversal__* of the given curves, which is computed by the function optimal_traversal() in the same file. Finally, the mean_of_curves() function returns *__the mean curve  of multiple given curves__*. <br><br>
  __Notable implementation decisions__:
  * The algorithm to compute the optimal traversal returns it in reverse (because vectors do not and should not have a way to push an element to the front). We decided not to reverse the returned traversal in order to not add additional computational weight to the program.
  * εδώ για το mean_of_curves

## Final thoughts




Οι δομές μας:
Για την αναπαράσταση των χρονοσειρών των δεδομένων μας ως διανύσματα χρησιμοποιούμε την κλάση Item (την οποία χρησιμοποιήσαμε και στο Project 1) η οποία χαρακτηρίζεται απο ένα id string και ένα vector απο doubles.
Για την αναπαράσταση των χρονοσειρών ως πολυγωνικές καμπύλες χρησιμοποιούμε την κλάση Curve2d. Η κλάση αυτή χαρακτηρίζεται μεταξύ άλλων απο ένα id string και έναν vector απο Point2d αντικείμενα. Τα point2d αντικείμενα αναπαριστούν ένα σημείο 2 διαστάσεων.
Για την τελική αποθήκευση των καμπυλών στα buckets χρησιμοποιούμε αντικείμενα κλάσης Association. Η κλάση αναπαριστά τη συσχέτιση μεταξύ
α) Της αρχικής καμπύλης
β) Της καμπύλης πλέγματος (grid curve) που λειτουργεί ώς το κλειδί του locality sensitive hashing
γ) Του x_vector που αναπαριστά την καμπύλη πλέγματος ως διάνυσμα και λειτουργεί ως κλειδί για το hashing της αποθήκευσης.

Η μέση καμπύλη ν καμπυλών υπολογίζεται χρησιμοποιώντας τη συνάρτηση mean_of_curves που προσομοιώνει έναν αλγόριθμο δυαδικού δέντρου ως εξής:
Η συνάρτηση παίρνει ως όρισμα ένα vector καμπυλών κλάσης Curve2d
Στο πρώτο στάδιο, ξεκινώντας απο τις 2 τελευταίες καμπύλες του vector, υπολογίζουμε το vector απο Point2d σημεία που αναπαριστά τη mean των δύο καμπυλών. Η mean του ζευγαριού γίνεται push σε έναν vector και ο δείκτης μειώνεται κατα 2. Αυτό γίνεται επαναληπτικά μέχρι 
α) είτε ο δείκτης να φτάσει στο 0:
Άρα έχουμε υπολογίσει τις μέσες καμπύλες κάθε διαδοχικού ζευγαριού.
β) ο δείκτης να φτάσει στο 1:
Άρα έχει μείνει μια καμπλύλη που δεν είχε ζευγάρι γιατί το συνολικό πλήθος των καμπυλών ήταν περιττό. Σε αυτήν την περίπτωση η καμπύλη γίνεται push στην δομή με τις mean καμπύλες των ζευγαριών.
Αφού έχουμε "φορτώσει" τα δεδομένα σε μια δομή διαφορετική της αρχικής η διαδικασία αλλάζει ως εξής:
1) οι mean καμπύλες των διαδοχικών ζευγαριών των καμπυλών γίνονται push στο τέλος της ίδιας δομής και όχι σε διαφορετική
2) ο δείκτης μειώνεται κατα 2 ώστε να μην πιάσει τις νεοεισαχθέντες mean καμπύλες στο βήμα/πέρασμα κατα το οποίο μπήκαν στη δομή
3) το ζευγάρι των καμπυλών του οποίου η mean καμπύλη υπολογίστηκε διαγράφεται απο τη δομη
Η παραπάνω διαδικασία εφαρμόζεται επαναληπτικά εως ότου μείνει μόνο 1 καμπύλη στη δομή η οποία θα αποτελέι και τη μέση καμπύλη όλων των καμπυλών στις οποίες εφαρμόστηκε ο αλγόριθμος. 

TimeSeries-ANN/ContinuousFrechet/cont_Frechet.hpp
Αρχείο επικεφαλίδας που περιέχει συναρτήσεις τύπου "black box" που 
α) μετατρέπουν τα αντικείμενα των κλάσεών μας σε αντικείμενα των κλάσεων που χρησιμοποιούνται στον φάκελο Fred που μας υποδείχθηκε να ενσωματώσουμε χωρίς αλλαγές
β) χρησιμοποιούν τις συναρτήσεις που υπάρχουν υλοποιημένες στον φάκελο Fred για τον υπολογισμό της Continuous Frechet Distance
Επίσης στο αρχείο περιέχεται συνάρτησει που βρίσκει με εξαντλητική αναζήτηση την καμπύλη (απο ένα συνόλο καμπυλών) με την μικρότερη Continuous Frechet distance απο μια καμπύλη.

TimeSeries-ANN/ContinuousFrechet/contF_LSH.hpp
Αρχείο επικεφαλίδας που αφορά το LSH hashing πολυγωνικών καμπυλών στην ευθεία R και την αναζήτηση κοντινότερου γείτον χρησιμοποιώντας τη μετρική Continuous Frechet distance

Evaluation
Συμπληρωματικά των μετρικών που ζητούνται απο την εκφώνηση για το evaluation των αλγορίθμων χρησιμοποιούμε και δυο δικές μας μετρικές που εκτυπώνονται στο std out
α) tApproximateAverage / tTrueAverage
β) average (distanceApproximate / distanceTrue)

Πορίσματα - Παρατηρήσεις - Ανάλυση
Ai) Κάθε χρονοσειρά αναπαρίσταται ως διάνυσμα - Χρησιμοποιούμε τους αλγορίθμους vector LSH και Hypercube της 1ης εργασίας. 

Οι αλγόριθμοι της 1ης εργασίας όπως αναμένεται δεν είναι ιδανικοί για τα δεδομένα των χρονοσειρών καθώς προορίζονται για πιο "απλά" δεδομένα. Βέβαια δεν είναι τελείως μη αξιποιήσιμοι. Μπορούμε σε λίγο μικρότερο χρόνο σε σχέση με την εξαντλητική αναζήτηση να βρούμε σχετικά εύστοχα αποτελέσματα σε κοντινούς γείτονες.

Aiii) Κάθε χρονοσειρά αναπαρίσταται ως πολυγωνική καμπύλη στην ευθεία R - LSH για Continuous Frechet

Η χρήση των ήδη υλοποιημένων συναρτήσεων που μας δόθηκαν μέσω "black box" συναρτήσεων, σε συνδυασμό με τη φύση της μετρικής καθιστούν τον υπολογισμό της Continuous Frechet πολύ αργό σε σχέση με τη discrete frechet. Ακόμα και με αρκετό φιλτράρισμα των καμπυλών για μείωση της πολυπλοκότητας, ο αλγόριθμος καθυστερεί πολύ λόγω της αργής πράξης του υπολογισμού. Παράλληλα η καθυστέρηση στο brute force καθιστά τον έλεγχο και την αξιολόγηση των αποτελεσμάτων χρονοβόρα υπόθεση. Ωστόσο κάναμε δοκιμές και διαπιστώσαμε ότι απο άποψη ακρίβειας ο αλγόριθμος είναι άριστος και μπορεί να ανταποκριθεί στα δεδομένα των χρονοσειρών.

Aii) Η αναπαράσταση των χρονοσειρών ως πολυγωνικές καμπύλες στον R^2, η χρήση της Discrete Frechet για την αναζήτηση του κοντινότερου γείτονα και ο αντίστοιχος LSH αλγόριθμος.

Αναμφίβολα έδειξαν την καλύτερη επίδοση για τα δεδομένα των χρονοσειρών. Μπορούμε σε πολύ μικρότερο χρόνο σε σχέση με την εξαντλητική αναζήτηση να βρούμε αρκετά εύστοχα και αξιόπιστα αποτελέσματα.

Querying Trick (για Frechet)
Υλοποιήθηκε σύμφωνα με τις διαφάνειες ώστε τα στοιχεία στο bucket που μας έστειλε το hashing του Query να ελέγχονται πρώτα μια φορά ωστε άν ενα στοιχείο δεν έχει identical grid curve με το στοιχείο της αναζήτησης να μην θεωρείται υποψήφιο ώστς να αποφεύγεται ο υπολογισμός της frechet distance του.
Ωστόσο δεν διαπιστώσαμε κάποια σημαντική βελτίωση, τουλάχιστον για τα δεδομένα μας, ειδικά στη discrete frechet περίπτωση όπου ο υπολογισμός της είναι αρκετά γρήγορος. Άρα εκεί την έχουμε απενεργοποιημένη ως default. Εκτιμούμε ότι στην continuous frechet περίπτωση μπορεί να έχει σημαντικό αντίκτυπο αλλά η δοκιμή και το evaluation ήταν δύσκολο ώστε να έχουμε concrete συμπέρασμα.

Delta Tuning
Η παράμετρος delta όταν δεν δίνεται απο τον χρήστη αυτορυθμίζεται με βάση τα δεδομένα ωστε να προσεγγίζει τη μέση απόσταση μεταξύ των διαδοχικών κορυφών των καμπυλών του input dataset. Η συνάρτηση για ένα αριθμό επαναλήψεων επιλέγει τυχαία ζευγάρια καμπυλών και αθροίζει τη μέση απόσταση των κορυφών και βρίσκει τη μέση τιμή της για όλα τα ζευγάρια. Η υλοποίηση βρίσκεται στο αρχείο curves.hpp

Padding
Προσθέτουμε συντεταγμένες με τιμή 10.000 (ένας μεγάλος αριθμός πολύ εκτός του πεδίου τιμών) σε έναν x_vector μέχρι να φτάσει μήκος ίσο της αρχικής καμπύλης στην ευθέια R ή διπλάσσιο της αρχικής 2d καμπύλης. Έτσι όλα τα x_vectors έχουν ίδιο μήκος ώστε να δουλέυει το LSH της 1ης εργασίας για αποθληκευση στα buckets.

Parameter Tuning (best parameters for each algorithm)

Vector_LSH L=1, k=1, tablesize=datasetsize/2 w = average L2 distance between dataset curves (as vectors)

Hypercube k=1. M=30, probes=1, w = average L2 distance between dataset curves (as vectors)

discrete_LSH L = 6, delta = avg dist between curve vertices, 
tablesize = dataset_size/8, threshold = dataset_size/4, querying trick = false

continuous_LSH L = 1, delta = avg dist between curve vertices, 
tablesize = dataset_size/4, threshold = dataset_size/4, querying trick = true

Ενδεικτικές Εκτελέσεις **************δε μου τρεχουν ετσι***************
$./bin/search –i datasets/nasd_input.csv –q datasets/nasd_query.csv –k 1 -L 1 -ο outputs/output.txt -algorithm LSH

$./bin/search –i datasets/nasd_input.csv –q datasets/nasd_query.csv –k 1 -M 30 -probes 1 -ο outputs/output.txt -algorithm Hypercube

$./bin/search –i datasets/nasd_input.csv –q datasets/nasd_query.csv -L 6 -ο outputs/output.txt -algorithm Frechet -metric discrete

$./bin/search –i datasets/nasd_input.csv –q datasets/nasd_query.csv -ο outputs/output.txt -algorithm Frechet -metric continuous





