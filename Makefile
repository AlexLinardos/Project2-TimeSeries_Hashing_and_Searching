CC=g++
CFLAGS=-O3

search: LSH_main.o
	$(CC) LSH_main.o -o bin/search $(CFLAGS)

LSH_main.o: ./TS-NN/Basic/LSH_main.cpp ./TS-NN/Basic/LSH.hpp ./TS-NN/Basic/NN_interface.hpp
	$(CC) -c ./TS-NN/Basic/LSH_main.cpp $(CFLAGS)

clean:
	rm bin/*
	rm *.o