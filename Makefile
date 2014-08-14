#
# Makefile for STLice
#

BUILD_DIR = build
SRC_DIR = .
LIBS_DIR = .

BUILD_TYPE = RELEASE

VERSION ?= DEV
CXX ?= gcc
#CFLAGS += -c -Wall -Wextra -Wold-style-cast -Woverloaded-virtual -std=c++11 -DVERSION=\"$(VERSION)\" -isystem libs
CFLAGS += -c -Wall -Wextra -Wold-style-cast -Woverloaded-virtual

ifeq ($(BUILD_TYPE),DEBUG)
	CFLAGS+=-ggdb -Og -g
endif
ifeq ($(BUILD_TYPE),PROFILE)
	CFLAGS+= -pg
endif
ifeq ($(BUILD_TYPE),RELEASE)
	CFLAGS+= -O3 -fomit-frame-pointer
endif

LDFLAGS += -Lbuild/ 	

SOURCES_RAW = main.c loadSTL.c 
SOURCES = $(addprefix $(SRC_DIR)/,$(SOURCES_RAW))

OBJECTS_RAW = $(SOURCES_RAW:.c=.o)
OBJECTS = $(addprefix $(BUILD_DIR)/,$(OBJECTS_RAW))

DIRS = $(sort $(dir $(OBJECTS)))

EXECUTABLE = $(BUILD_DIR)/STLice

ifeq ($(OS),Windows_NT)
	#For windows make it large address aware, which allows the process to use more then 2GB of memory.
	EXECUTABLE := $(EXECUTABLE).exe
	CFLAGS += -march=pentium4 -flto
	LDFLAGS += -Wl,--large-address-aware -lm -lwsock32 -flto
	MKDIR_PREFIX = mkdir -p
else
	MKDIR_PREFIX = mkdir -p
	UNAME := $(shell uname)
	ifeq ($(UNAME), Linux)
		OPEN_HTML=firefox
		ifeq ($(BUILD_TYPE),DEBUG)
			LDFLAGS += --static
		else
			CFLAGS += -flto
			LDFLAGS += --static -flto
		endif
	endif
	ifeq ($(UNAME), Darwin)
		OPEN_HTML=open
		#For MacOS force to build
		CFLAGS += -force_cpusubtype_ALL -mmacosx-version-min=10.6 -arch x86_64 -arch i386
		LDFLAGS += -force_cpusubtype_ALL -mmacosx-version-min=10.6 -arch x86_64 -arch i386
	endif
endif

all: $(DIRS) $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

$(DIRS):
	-@$(MKDIR_PREFIX) $(DIRS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CXX) $(CFLAGS) $< -o $@

## clean stuff
clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

help:
	@cat Makefile |grep \#\#| grep \: |cut -d\# -f3
