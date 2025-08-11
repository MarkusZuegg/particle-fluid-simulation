CC = g++

CFLAGS = -Wall -pedantic -Wextra -I"${HDF5_INSTALL_PATH}/include"

LDFLAGS = -L"${HDF5_INSTALL_PATH}/lib" -lhdf5_cpp -lm

DEBUG = -g

.DEFAULT_GOAL := all

.PHONY: clean

debug: CFLAGS += $(DEBUG)
debug: clean $(TARGET)

all: sph

sph: sph.cpp
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f sph *.o
