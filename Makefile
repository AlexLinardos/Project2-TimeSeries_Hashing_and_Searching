CC=g++
CFLAGS=-O3
CXXFLAGS = -march=native -Ofast -static-libgcc -static-libstdc++ -std=c++14 -fpermissive -fPIC -ffast-math -fno-trapping-math -ftree-vectorize

all: search clean1

search: search_main.o config.o curve.o frechet.o interval.o point.o simplification.o
	$(CC) search_main.o config.o curve.o frechet.o interval.o point.o simplification.o -o bin/search $(CFLAGS)

search_main.o: ./TS-NN/search_main.cpp ./TS-NN/Basic/LSH.hpp ./TS-NN/Basic/HC.hpp ./TS-NN/NN_interface.hpp
	$(CC) -c ./TS-NN/search_main.cpp $(CFLAGS) $(CXXFLAGS)

config.o: ./TS-NN/ContinuousFrechet/Fred/config.cpp ./TS-NN/ContinuousFrechet/Fred/config.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/config.cpp $(CFLAGS) $(CXXFLAGS)

curve.o: ./TS-NN/ContinuousFrechet/Fred/curve.cpp ./TS-NN/ContinuousFrechet/Fred/curve.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/curve.cpp $(CFLAGS) $(CXXFLAGS)

frechet.o: ./TS-NN/ContinuousFrechet/Fred/frechet.cpp ./TS-NN/ContinuousFrechet/Fred/frechet.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/frechet.cpp $(CFLAGS) $(CXXFLAGS)

interval.o: ./TS-NN/ContinuousFrechet/Fred/interval.cpp ./TS-NN/ContinuousFrechet/Fred/interval.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/interval.cpp $(CFLAGS) $(CXXFLAGS)

point.o: ./TS-NN/ContinuousFrechet/Fred/point.cpp ./TS-NN/ContinuousFrechet/Fred/point.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/point.cpp $(CFLAGS) $(CXXFLAGS)

simplification.o: ./TS-NN/ContinuousFrechet/Fred/simplification.cpp ./TS-NN/ContinuousFrechet/Fred/simplification.hpp ./TS-NN/ContinuousFrechet/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFrechet/Fred/simplification.cpp $(CFLAGS) $(CXXFLAGS)

clean:
	rm bin/*

clean1:
	rm *.o