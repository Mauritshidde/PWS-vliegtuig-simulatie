objects = main.o simulation.o
CXX = g++

all: main

main: ${objects}
	${CXX} -o main ${objects} -lraylib

simulation.o: simulation.cpp simulation.h ModelLoader.o
	${CXX} -c simulation.cpp -lraylib

main.o: main.cpp
	${CXX} -c main.cpp

ModelLoader.o: ModelLoader.cpp ModelLoader.h
	${CXX} -c ModelLoader.cpp
clean:
	rm *.o main