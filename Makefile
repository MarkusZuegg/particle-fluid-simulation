# Directories
SRCDIR = src
LIBDIR = lib
BINDIR = bin

# Compilers
CC = gcc
CXX = g++
NVCC = nvcc

CFLAGS = -std=c11
LIBCFLAGS = -std=c11
CXXFLAGS = -Wall -pedantic -Wextra -std=c++2a -I$(LIBDIR) -L$(LIBDIR)

DEBUG = -g
OPT = -O2

CUDA_INC = /opt/local/stow/cuda-11.1/include
CUDA_LIB = /opt/local/stow/cuda-11.1/lib64
GPU_ARCH = sm_35
NVCCFLAGS = -std=c++17 -arch=$(GPU_ARCH) -I$(CUDA_INC) -I$(LIBDIR)

LDFLAGS = -L$(LIBDIR) -L$(CUDA_LIB) -lcudart

TARGET_CPU = $(SRCDIR)/sph_cpu
TARGET_GPU = $(SRCDIR)/sph_gpu
TARGET_GPU_SOA = $(SRCDIR)/sph_gpu_SoA

SRC_C = $(wildcard $(LIBDIR)/*.c)
OBJS_C = $(patsubst $(LIBDIR)/%.c,$(BINDIR)/%.o,$(SRC_C))

OBJS_CPU = $(BINDIR)/sph_cpu.o $(OBJS_C)
OBJS_GPU = $(BINDIR)/sph_gpu.o $(OBJS_C)
OBJS_GPU_SOA = $(BINDIR)/sph_gpu_SoA.o $(OBJS_C)

# Default goal
.DEFAULT_GOAL := all
.PHONY: all cpu gpu clean debug opti gpu_SoA

all: $(TARGET_CPU) $(TARGET_GPU)

debug: CFLAGS += $(DEBUG)
debug: NVCCFLAGS += $(DEBUG)
debug: CXXFLAGS += $(DEBUG)
debug: clean $(TARGET)

opti: CFLAGS += $(OPT)
opti: NVCCFLAGS += $(OPT)
opti: LIBCFLAGS += $(OPT)
opti: CXXFLAGS += $(OPT)
opti: clean $(TARGET)

cpu: $(TARGET_CPU)
gpu: $(TARGET_GPU)
gpu_SoA: $(TARGET_GPU_SOA)

$(TARGET_CPU): $(OBJS_CPU)
	$(CXX) $(CXXFLAGS) $(OBJS_CPU) -o $@

$(TARGET_GPU): $(OBJS_GPU)
	$(NVCC) $(NVCCFLAGS) $(OBJS_GPU) -o $@ $(LDFLAGS)

$(TARGET_GPU_SOA): $(OBJS_GPU_SOA)
	$(NVCC) $(NVCCFLAGS) $(OBJS_GPU_SOA) -o $@ $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/sph_cpu.o: $(SRCDIR)/sph.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BINDIR)/sph_gpu.o: $(SRCDIR)/sph-cuda.cu | $(BINDIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BINDIR)/sph_gpu_SoA.o: $(SRCDIR)/sph-cuda-SoA.cu | $(BINDIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(LIBDIR)/%.c | $(BINDIR)
	$(CC) $(LIBCFLAGS) -c $< -o $@

clean:
	rm -rf $(BINDIR) $(TARGET_CPU) $(TARGET_GPU)
