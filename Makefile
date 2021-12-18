CC=g++
CFLAGS=-O3
CXXFLAGS = -march=native -Ofast -static-libgcc -static-libstdc++ -std=c++14 -fpermissive -fPIC -ffast-math -fno-trapping-math -ftree-vectorize

all: search cluster

search: final_search clean1

final_search: search_main.o config.o curve.o frechet.o interval.o point.o simplification.o
	$(CC) search_main.o config.o curve.o frechet.o interval.o point.o simplification.o -o bin/search $(CFLAGS)

search_main.o: ./src/search_main.cpp ./TimeSeries-ANN/L2/LSH.hpp ./TimeSeries-ANN/L2/HC.hpp ./ui/NN_interface.hpp
	$(CC) -c ./src/search_main.cpp $(CFLAGS) $(CXXFLAGS)

config.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/config.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/config.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/config.cpp $(CFLAGS) $(CXXFLAGS)

curve.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/curve.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/curve.hpp ./TimeSeries-ANN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/curve.cpp $(CFLAGS) $(CXXFLAGS)

frechet.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/frechet.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/frechet.hpp ./TimeSeries-ANN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/frechet.cpp $(CFLAGS) $(CXXFLAGS)

interval.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/interval.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/interval.hpp ./TimeSeries-ANN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/interval.cpp $(CFLAGS) $(CXXFLAGS)

point.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/point.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/point.hpp ./TimeSeries-ANN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/point.cpp $(CFLAGS) $(CXXFLAGS)

simplification.o: ./TimeSeries-ANN/ContinuousFrechet/Fred/simplification.cpp ./TimeSeries-ANN/ContinuousFrechet/Fred/simplification.hpp ./TimeSeries-ANN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TimeSeries-ANN/ContinuousFrechet/Fred/simplification.cpp $(CFLAGS) $(CXXFLAGS)

cluster: final_cluster clean1

final_cluster: cluster_main.o
	$(CC) cluster_main.o -o bin/cluster $(CFLAGS)

cluster_main.o: ./src/cluster_main.cpp ./ui/Clustering_interface.hpp
	$(CC) -c ./src/cluster_main.cpp $(CFLAGS) $(CXXFLAGS)

clean:
	rm bin/*
	rm outputs/*

clean1:
	rm *.o
