# ============================================================================
#  Makefile for the UParT b-tag data/MC plotting code
#
#  On this machine ROOT is installed as a snap, which brings its own compiler.
#  Run all make commands inside the snap environment:
#
#      root-framework.run make          # compile
#      root-framework.run make run      # full run
#      root-framework.run make clean    # remove compiled files
#
#  (On a machine with a normal ROOT installation, plain `make` works.)
# ============================================================================

# Compiler and flags. root-config tells us where ROOT's headers and libraries are.
CXX      := g++
CXXFLAGS := -O2 -g -Wall -Wextra $(shell root-config --cflags) -Iinclude
LDFLAGS  := $(shell root-config --ldflags)
LDLIBS   := $(shell root-config --libs)

# Every .cc file in src/ is compiled into build/<name>.o and then linked
# together into one program, bin/makePlots.
# SRCS   := $(wildcard src/*.cc)
# OBJS   := $(SRCS:src/%.cc=build/%.o)
# Nur die spezifischen Quelldateien für unser Programm definieren
SRCS   := src/Plotting.cc src/make2DPlots.cc
OBJS   := build/Plotting.o build/make2DPlots.o
DEPS   := $(OBJS:.o=.d)
TARGET := bin/makePlots

# Fail early with a helpful message if the compiler is not found
ifeq ($(shell command -v $(CXX) 2>/dev/null),)
$(error Compiler '$(CXX)' not found. With the ROOT snap, use 'root-framework.run make' instead of 'make')
endif

.PHONY: all run clean

all: $(TARGET)

# Link
$(TARGET): $(OBJS) | bin
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compile. -MMD -MP also writes build/<name>.d, which lists the headers each
# file includes, so that changing a header recompiles the right files.
build/%.o: src/%.cc | build
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

bin build:
	mkdir -p $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf build bin

-include $(DEPS)
