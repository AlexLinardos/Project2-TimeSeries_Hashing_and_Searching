# Project2-TimeSeries_Hashing_and_Searching

## Runs
./search -> Τρέχει με τις παραμέτρους καρφωτές από το NN_interface (γραμμές 222-231). Το έκανα για να το τρέχω κατευθείαν από vs code με το κουμπί <br>

./search -i ../datasets/nasd_input.csv -q ../datasets/nasd_query.csv -o output.txt -algorithm LSH -delta 0.0 -> Αυτές είναι οι υποχρεωτικές παράμετροι. Όλες οι άλλες αν δεν δοθούν θα αρχικοποιηθούν με default τιμή βάσει εκφώνησης. <b>Αν -algorithm Frechet τότε είναι υποχρεωτική και η -metric.</b> Δηλαδή:<br>

./search -i ../datasets/nasd_input.csv -q ../datasets/nasd_query.csv -o output.txt -algorithm Frechet -metric discrete -delta 0.0 <br>

<b>Note: Κανονικά ούτε η -delta είναι υποχρεωτική, αλλά απ ότι κατάλαβα θα πρέπει να υπολογίσουμε εμείς (ή να υπολογίζει το πρόγραμμα καλύτερα) την τιμή της βάσει του dataset οπότε προς το παρόν την έχω αφήσει.</b>