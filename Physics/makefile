objects = matrix.o navierstokes.o
CXX = g++

all: navierstokes

navierstokes: ${objects}
	${CXX} -o navierstokes ${objects} -I /usr/include/python3.10/ -lpython3.10

# matplotlib.o: matplotlib.h
# 	${CXX} -c matplotlib.h

matrix.o: matrix.cpp matrix.h
	${CXX} -c matrix.cpp

navierstokes.o: navierstokes.cpp navierstokes.h
	${CXX} -c navierstokes.cpp -I /usr/include/python3.10/ -lpython3.10

clean:
	rm *.o navierstokes