# objects = main.o simulation

# all: main

# main: ${objects}
# 	${CXX} -o main ${objects} -lraylib

# ModelLoader.o: ModelLoader.cpp ModelLoader.h
# 	${CXX} -c ModelLoader.cpp

# simulation: simulation.o ModelLoader.o
# 	${CXX} -c simulation.o ModelLoader.o

# simulation.o: simulation.cpp simulation.h
# 	${CXX} -c simulation.cpp
# main.o: main.cpp
# 	${CXX} -c main.cpp

# clean:
# 	rm *.o main
TARGET ?= a.out
SRC_DIRS ?= ./src
CC = g++

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d) /usr/include/python3.11 
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@ $(LOADLIBES) $(LDLIBS) -lraylib

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)