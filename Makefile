CC=g++
CFLAGS=-O3
CXXFLAGS = -march=native -Ofast -static-libgcc -static-libstdc++ -std=c++14 -fpermissive -fPIC -ffast-math -fno-trapping-math -ftree-vectorize

all: search clean1

search: search_main.o
	$(CC) search_main.o config.o curve.o frechet.o interval.o point.o simplification.o -o bin/search $(CFLAGS)

search_main.o: ./TS-NN/search_main.cpp ./TS-NN/Basic/LSH.hpp ./TS-NN/Basic/HC.hpp ./TS-NN/NN_interface.hpp
	$(CC) -c ./TS-NN/search_main.cpp $(CFLAGS) $(CXXFLAGS)

config.o: ./TS-NN/ContinuousFreche/Fred/config.cpp ./TS-NN/ContinuousFreche/Fred/config.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/config.cpp $(CXXFLAGS)

curve.o: ./TS-NN/ContinuousFreche/Fred/curve.cpp ./TS-NN/ContinuousFreche/Fred/curve.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/curve.cpp $(CXXFLAGS)

frechet.o: ./TS-NN/ContinuousFreche/Fred/frechet.cpp ./TS-NN/ContinuousFreche/Fred/frechet.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/frechet.cpp $(CXXFLAGS)

interval.o: ./TS-NN/ContinuousFreche/Fred/interval.cpp ./TS-NN/ContinuousFreche/Fred/interval.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/interval.cpp $(CXXFLAGS)

point.o: ./TS-NN/ContinuousFreche/Fred/point.cpp ./TS-NN/ContinuousFreche/Fred/point.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/point.cpp $(CXXFLAGS)

simplification.o: ./TS-NN/ContinuousFreche/Fred/simplification.cpp ./TS-NN/ContinuousFreche/Fred/simplification.hpp ./TS-NN/ContinuousFreche/Fred/types.hpp
	$(CC) -c ./TS-NN/ContinuousFreche/Fred/simplification.cpp $(CXXFLAGS)

clean:
	rm bin/*

clean1:
	rm *.o