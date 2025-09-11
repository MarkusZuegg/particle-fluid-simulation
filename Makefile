# Directories
SRCDIR = src
LIBDIR = lib
BINDIR = bin

# Compilers
CC = gcc
CXX = g++
CFLAGS = -std=c11
LIBCFLAGS = -std=c11
CXXFLAGS = -Wall -pedantic -Wextra -std=c++2a -I$(LIBDIR) -L$(LIBDIR)
DEBUG = -g

# Target
TARGET = $(SRCDIR)/sph

# Sources
SRC_CPP = $(SRCDIR)/sph.cpp
SRC_C = $(wildcard $(LIBDIR)/*.c)

# Object files in bin/
OBJS = $(BINDIR)/sph.o $(patsubst $(LIBDIR)/%.c,$(BINDIR)/%.o,$(SRC_C))

# Default goal
.DEFAULT_GOAL := all
.PHONY: clean debug all

all: $(TARGET)

debug: CFLAGS += $(DEBUG)
debug: CXXFLAGS += $(DEBUG)
debug: clean $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Ensure bin directory exists
$(BINDIR):
	mkdir -p $(BINDIR)

# Compile C++ source
$(BINDIR)/sph.o: $(SRCDIR)/sph.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile C library sources
$(BINDIR)/%.o: $(LIBDIR)/%.c | $(BINDIR)
	$(CC) $(LIBCFLAGS) -c $< -o $@

clean:
	rm -rf $(BINDIR) $(TARGET)
