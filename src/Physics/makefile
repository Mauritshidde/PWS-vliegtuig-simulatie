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

	
debug:
	@echo " test$(python_version_full)"	
	@echo "test $(PYTHON_VERSION)"
	@echo "test"
	@echo "test"
	@echo "test"
	@echo "/usr/include/$(python_version_major).$(python_version_minor)"
	@echo ${python_version_full}
	@echo $(python_version_major)
	@echo ${python_version_minor}
	@echo ${python_version_patch}
	@echo ${my_cmd}