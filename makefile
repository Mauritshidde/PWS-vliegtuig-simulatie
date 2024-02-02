TARGET ?= a.out
SRC_DIRS ?= ./src
CC = nvcc
NVCCLFLAGS := -std=c++11 -cudart=shared -rdc=true

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.cu' -or -name '*.c' -or -name '*.s')
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
python_version_patch := $(word 3,${python_version_full})

INC_DIRS := $(shell find $(SRC_DIRS) -type d) /usr/include/python$(python_version_major).$(python_version_minor)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

# Specify the compilation rule for CUDA files
%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

%.o: %.cu
	$(CC) $(NVCCLFLAGS) $(CPPFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(NVCCLFLAGS) $(LDFLAGS) $(OBJS) -o $@ $(LOADLIBES) $(LDLIBS) -lraylib -lpython$(python_version_major).$(python_version_minor)

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)
