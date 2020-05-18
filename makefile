all: bin/omDBG

bin/omDBG:
	mkdir bin
	g++ -w -std=c++11 src/*.cpp -o bin/rmapper
	g++ -std=c++11 graph_traverse.cpp -o traverse
	

clean: 
	rm -rf bin