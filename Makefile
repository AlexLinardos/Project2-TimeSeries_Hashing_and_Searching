CC=g++
CFLAGS=-O3

all: search clean1

search: search_main.o
	$(CC) search_main.o -o bin/search $(CFLAGS)

search_main.o: ./TS-NN/search_main.cpp ./TS-NN/Basic/LSH.hpp ./TS-NN/Basic/HC.hpp ./TS-NN/NN_interface.hpp
	$(CC) -c ./TS-NN/search_main.cpp $(CFLAGS)

clean:
	rm bin/*

clean1:
	rm *.o