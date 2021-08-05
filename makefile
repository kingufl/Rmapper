all: rmapper

rmapper:
	g++ -w -std=c++11 src/*.cpp -o rmapper-build
	g++ -std=c++11 graph_traverse.cpp -o rmapper-traverse
	

clean: 
	rm -rf bin