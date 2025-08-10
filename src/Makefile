CC = g++

CFLAGS = -Wall -pedantic -std=gnu99 -Wextra -lhdf5

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
