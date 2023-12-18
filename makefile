objects = simulation.o
CXX = g++

all: main

main: ${objects}
	${CXX} -o main ${objects}

# matplotlib.o: matplotlib.h
# 	${CXX} -c matplotlib.h

simulation.o: simulation.cpp simulation.h
	${CXX} -c simulation.cpp

clean:
	rm *.o main