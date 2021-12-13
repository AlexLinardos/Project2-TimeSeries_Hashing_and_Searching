CC=g++
CFLAGS=-O3

search: search_main.o
	$(CC) search_main.o -o bin/search $(CFLAGS)

search_main.o: ./TS-NN/Basic/search_main.cpp ./TS-NN/Basic/LSH.hpp ./TS-NN/Basic/NN_interface.hpp
	$(CC) -c ./TS-NN/Basic/search_main.cpp $(CFLAGS)

clean:
	rm bin/*
	rm *.o