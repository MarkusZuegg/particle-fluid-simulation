CC = g++

CFLAGS = -Wall -pedantic -Wextra -std=c++17

LDFLAGS = -L"${HDF5_INSTALL_PATH}/lib" -lhdf5_cpp -lm

DEBUG = -g

.DEFAULT_GOAL := all

.PHONY: clean

debug: CFLAGS += $(DEBUG)
debug: clean $(TARGET)

all: sph

sph: sph.cpp
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f sph *.o
