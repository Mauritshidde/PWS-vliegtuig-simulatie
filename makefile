TARGET ?= a.out
SRC_DIRS ?= ./src
CC = g++

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.s')
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
python_version_patch := $(word 3,${python_version_full})

INC_DIRS := $(shell find $(SRC_DIRS) -type d) /usr/include/python$(python_version_major).$(python_version_minor)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@ $(LOADLIBES) $(LDLIBS) -lraylib -lpython3.10

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)