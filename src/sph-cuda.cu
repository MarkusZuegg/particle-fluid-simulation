// Created by Markus Zuegg
// Main cpp script to run sph simulation
// Main simulation loop will be housed here
#include <bits/stdc++.h>
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include "cJSON.h"
#include "gsd.h"
#include <math.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/for_each.h>
#include <thrust/sort.h>
#include <thrust/fill.h>
#include <thrust/iterator/counting_iterator.h>

using std::vector;
using timer = std::chrono::high_resolution_clock;
using mircoseconds = std::chrono::microseconds;

// =================== Global Variables ========================================
const int numDimensions = 2;

// =========================== Structs =========================================
// Bounding box for simulation
// x min and max
// y min and max
typedef struct
{
    float xmax;
    float xmin;
    float ymax;
    float ymin;
    float zmax;
    float zmin;
} box;

// Holds all the relivant params needed for the simulation
typedef struct
{
    int numParticles;
    float timeStep;
    float gravity;
    float elasticity;
    int mass;
    float radius;
    float smoothingRadius;
    float targetDensity;
    float pressureMulti;
    float viscosity;
} SimulationParams;


// Struct to hold all varibales needed for a particle
// AoS implementation
typedef struct
{
    float x;
    float y;
    float vx;
    float vy;
    float density;
    int gridHash;
} Particle;

//! currently depreicated (maybe use later)
// Ensum to hold all version that can be run
enum Versions {
    SoA = 0, //n^2 approach
    SoA_LINKED = 1, // struct of arrays, with linkedlist grid
    AoS_LINKED = 2, // array of structs, with linkedlist grid
    AoS_SORTED = 3 // array of structs, with sorted grid
};

// ========================= Func Protos =======================================

__global__ void find_cell_start_kernal(const Particle* particleArray, int* cellStart, int numParticles);
__device__ inline float ploy6_kernal(float dist, float radius);
__device__ inline float spikey_kernal(float dist, float radius);
__device__ inline float spikey_kernal_derivative(float dist, float radius);
__device__ inline float viscosity_kernal(float dist, float radius);
__device__ float distance_from_two_points(float x1, float y1, float x2, float y2);
__device__ float density_to_pressure(float density, const SimulationParams &params);

// General
int get_x_index(int i);
int get_y_index(int i);

// Inits
void init_random_particle_pos_AoS(vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params);
// Adjusted for thrust vector
void init_grid_particle_pos_AoS(thrust::host_vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params);

// File saving
void save_to_json(const std::string &filename, const SimulationParams &params, const box &bounds);
void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<float> &positions, const vector<float> &velocities,
        const SimulationParams &params);

// ===================== Thrust Functors =======================

struct ApplyHalfGravityFunctor
{
    const SimulationParams params;
    __host__ ApplyHalfGravityFunctor(const SimulationParams& params) : params(params) {}
    __device__ void operator()(Particle& p) const {
        p.vy += params.gravity * params.timeStep * 0.5f;
    }
};

// Given a particle calculate its grid hash
// @param SimulationParams params
// @param Box bounds
// @param float cellSize
// @param int gridRows
// @param int gridCols
struct BuildGridHashFunctor
{
    const SimulationParams params;
    const box bounds;
    const float cellSize;
    const int gridRows;
    const int gridCols;

    __host__ BuildGridHashFunctor(const SimulationParams& params, const box& b, float cs, int gr, int gc)
        : params(params), bounds(b), cellSize(cs), gridRows(gr), gridCols(gc) {}

    __device__ void operator()(Particle& p) const {
        int cellX = (int)floorf((p.x - bounds.xmin) / cellSize);
        int cellY = (int)floorf((p.y - bounds.ymin) / cellSize);
        // clamp with out using std::clamp
        if (cellX < 0) cellX = 0;
        if (cellX >= gridCols) cellX = gridCols - 1;
        if (cellY < 0) cellY = 0;
        if (cellY >= gridRows) cellY = gridRows - 1;

        p.gridHash = cellX + (cellY * gridCols);
    }
};

struct ParticleHashComparator
{
    __host__ __device__ bool operator()(const Particle& a, const Particle& b) const {
        return a.gridHash < b.gridHash;
    }
};

// Calculate density of a particle
// @param SimulationParams params
// @param Box bounds
// @param float cellSize
// @param int gridRows
// @param int gridCols
// @param int* cellStart
// @param Particle* particleArray
struct CalcDensityFunctor
{
    const SimulationParams params;
    const box bounds;
    const float cellSize;
    const int gridRows;
    const int gridCols;
    const int* cellStart;
    Particle* particleArray;

    __host__ CalcDensityFunctor(const SimulationParams& params, const box& b, float cs, int gr, int gc,
        const int* cs_ptr, Particle* p_ptr)
        : params(params), bounds(b), cellSize(cs), gridCols(gc), gridRows(gr), cellStart(cs_ptr), 
        particleArray(p_ptr) {}

    __device__ void operator()(int i) const {
        float density = 0.0;
        float x1 = particleArray[i].x;
        float y1 = particleArray[i].y;

        int cellX = static_cast<int>(floorf((x1 - bounds.xmin) / cellSize));
        int cellY = static_cast<int>(floorf((y1 - bounds.ymin) / cellSize));
        // clamp with out using std::clamp
        if (cellX < 0) cellX = 0;
        if (cellX >= gridCols) cellX = gridCols - 1;
        if (cellY < 0) cellY = 0;
        if (cellY >= gridRows) cellY = gridRows - 1;
        
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int newCellX = cellX + dx;
                int newCellY = cellY + dy;
                if (newCellX < 0 || newCellX >= gridCols) continue;
                if (newCellY < 0 || newCellY >= gridRows) continue;

                int neighbourCell = newCellX + (newCellY * gridCols);
                int start_idx = cellStart[neighbourCell];

                // If cell is empty
                if (start_idx == -1) continue;

                // Loop untill particle's hash != current cell
                for (int j = start_idx; j < params.numParticles && particleArray[j].gridHash == neighbourCell; j++) {
                    float x2 = particleArray[j].x;
                    float y2 = particleArray[j].y;
                    float dist = distance_from_two_points(x1, y1, x2, y2);
                    float influence = ploy6_kernal(dist, params.smoothingRadius);
                    density += params.mass * influence;
                }
            }
        }
        particleArray[i].density = density;
    }
};

// Calculate pressure forces on a particle
// @param SimulationParams params
// @param Box bounds
// @param float cellSize
// @param int gridRows
// @param int gridCols
// @param int* cellStart
// @param Particle* particleArray
struct CalcPressureFunctor
{
    const SimulationParams params;
    const box bounds;
    const float cellSize;
    const int gridRows;
    const int gridCols;
    const int* cellStart;
    Particle* particleArray;

    __host__ CalcPressureFunctor(const SimulationParams& params, const box& b, float cs, int gr, int gc,
        const int* cs_ptr, Particle* p_ptr)
        : params(params), bounds(b), cellSize(cs), gridCols(gc), gridRows(gr), cellStart(cs_ptr), 
        particleArray(p_ptr) {}

    __device__ void operator()(int i) const {
        float force_x = 0.0;
        float force_y = 0.0;

        float x1 = particleArray[i].x;
        float y1 = particleArray[i].y;

        // Get current particle velocity
        float vx1 = particleArray[i].vx;
        float vy1 = particleArray[i].vy;

        float density1 = particleArray[i].density;
        if (density1 <= 0.0) return;
        float pressure1 = density_to_pressure(density1, params);

        // get current cell pos
        int cellX = (int)floorf((x1 - bounds.xmin) / cellSize);
        int cellY = (int)floorf((y1 - bounds.ymin) / cellSize);
        // clamp with out using std::clamp
        if (cellX < 0) cellX = 0;
        if (cellX >= gridCols) cellX = gridCols - 1;
        if (cellY < 0) cellY = 0;
        if (cellY >= gridRows) cellY = gridRows - 1;

        // loop through surounding cells
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int newCellX = cellX + dx;
                int newCellY = cellY + dy;
                if (newCellX < 0 || newCellX >= gridCols) continue;
                if (newCellY < 0 || newCellY >= gridRows) continue;

                int neighbourCell = newCellX + (newCellY * gridCols);
                int start_idx = cellStart[neighbourCell];

                // If cell is empty
                if (start_idx == -1) continue;

                // Loop untill particle's hash != current cell
                for (int j = start_idx; j < params.numParticles && particleArray[j].gridHash == neighbourCell; j++) {
                    if (&particleArray[i] == &particleArray[j]) continue;

                    // get other particle postion
                    float x2 = particleArray[j].x;
                    float y2 = particleArray[j].y;

                    // Get distance/ magnitiude
                    float dist = distance_from_two_points(x1, y1, x2, y2);
                    if (dist == 0) continue; // Either particle is same, or ontop of eachother

                    // Get unit vectors
                    float dirX = (x1 - x2) / dist;
                    float dirY = (y1 - y2) / dist;

                    // Get the slope of the kernal at the distance
                    float slope = spikey_kernal_derivative(dist, params.smoothingRadius); // Grab the gradient of pressure change
                    // Get density and pressure for 2nd particle
                    float density2 = particleArray[j].density;
                    float pressure2 = density_to_pressure(density2, params);

                    float gradX = dirX * slope;
                    float gradY = dirY * slope;
                    // Symetric pressure calculation
                    float coeff = -params.mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

                    force_x += coeff * gradX;
                    force_y += coeff * gradY;

                    // do viscosity things
                    float influence = viscosity_kernal(dist, params.smoothingRadius);

                    float vx2 = particleArray[j].vx;
                    float vy2 = particleArray[j].vy;
                    
                    float vx_diff = vx2 - vx1;
                    float vy_diff = vy2 - vy1;

                    force_x += params.viscosity * vx_diff * influence;
                    force_y += params.viscosity * vy_diff * influence;
                }
            }
        }
        particleArray[i].vx += (force_x / density1) * params.timeStep;
        particleArray[i].vy += (force_y / density1) * params.timeStep;
    }
};

// Bounces particles back into box if needed
// @param SimulationParams params
// @param Box bounds 
struct BounceParticleFunctor
{
    const SimulationParams params;
    const box bounds;

    __host__ BounceParticleFunctor(const SimulationParams& params, const box& b) : params(params), bounds(b) {}
    
    __device__ void operator()(Particle& p) const {
        // check walls in x
        if (p.x - params.radius < bounds.xmin) {
            p.x = bounds.xmin + params.radius;
            p.vx *= -1 * params.elasticity;
        }
        if (p.x + params.radius > bounds.xmax) {
            p.x = bounds.xmax - params.radius;
            p.vx *= -1 * params.elasticity;
        }
        // Check walls in y
        if (p.y - params.radius < bounds.ymin) {
            p.y = bounds.ymin + params.radius;
            p.vy *= -1 * params.elasticity;
        }
        if (p.y + params.radius > bounds.ymax) {
            p.y = bounds.ymax - params.radius;
            p.vy *= -1 * params.elasticity;
        }
    }
};

// Updates the particles postion due to velocity
// @param SimulationParams params
struct UpdateParticleFunctor
{
    const SimulationParams params;

    __host__ UpdateParticleFunctor(const SimulationParams& params) : params(params) {}

    __device__ void operator()(Particle& p) const {
        p.x += p.vx * params.timeStep;
        p.y += p.vy * params.timeStep;
    }
};

// ========================== MAIN ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg which is jobID

    // Total time
    mircoseconds totalTime{0};

    // ============= Simulation params ==================
    const int maxNumSteps = 8000;
    const int saveFrames = 10; // How many frames to save. != 0
    const int numParticles = 1225; // 35 grid
    const float timeStep = 0.001;

    // gravity!!! HERE
    const float gravity = -9.81;
    const float elasticity = 0.75;

    const int mass = 1; // Particle mass
    const float radius = 0.1; // for particles

    const float smoothingRadius = 0.3;
    const float targetDensity = 20;
    const float pressureMulti = 1200;
    const float viscosity = 0.1;

    const box bounds = {12, 0, 10, 0, 0, 0}; //xmax, xmin, ymax, ymin, zmax, zmin

    // Change the line below to switch between each version
    int runType = AoS_SORTED;

    // Set up struct to hold params
    const SimulationParams params = {
        numParticles,
        timeStep,
        gravity,
        elasticity,
        mass,
        radius,
        smoothingRadius,
        targetDensity,
        pressureMulti,
        viscosity
    };

    // 2 array linked list setup
    float cellSize = smoothingRadius;
    int gridRows = std::max(1, (int)ceil((bounds.ymax - bounds.ymin) / cellSize));
    int gridCols = std::max(1, (int)ceil((bounds.xmax - bounds.xmin) / cellSize));
    int numCells = gridRows * gridCols;
    vector<int> cellStart;

    // ========================== File Names ==============================
    std::string fileNameCsv = "logs/x" + std::string(argv[0]) + ".csv";
    std::string fileNameJson = "logs/x" + std::string(argv[0]) + ".json";
    std::string fileNameGsd = "logs/x" + std::string(argv[0]) + ".gsd";

    // =========================== GSD ====================================
    gsd_handle handle;
    int ret = gsd_create_and_open(&handle, fileNameGsd.c_str(), "SPH simulation",
            "hoomd", gsd_make_version(1, 0), GSD_OPEN_READWRITE, 0);
    if (ret != 0) {
        std::cerr << "Error opening GSD file" << std::endl;
        return -1;
    }

    // ===================== Simulation loop ==================================
    thrust::host_vector<Particle> h_particleArray(numParticles);
    init_grid_particle_pos_AoS(h_particleArray, bounds, params);

    thrust::device_vector<Particle> d_particleArray = h_particleArray;
    thrust::device_vector<int> d_cellStart(numCells);

    for (int step = 0; step < maxNumSteps; step++) {
        // Time execution
        auto start = timer::now();
    
        // Apply half vy step as per leapfrog
        thrust::for_each(d_particleArray.begin(), d_particleArray.end(),
            ApplyHalfGravityFunctor(params));

        // Update position
        thrust::for_each(d_particleArray.begin(), d_particleArray.end(),
            UpdateParticleFunctor(params));

        // Build grid using thurt for each
        thrust::for_each(d_particleArray.begin(), d_particleArray.end(),
            BuildGridHashFunctor(params, bounds, cellSize, gridRows, gridCols));

        // Sort using thrust :)
        thrust::sort(d_particleArray.begin(), d_particleArray.end(), ParticleHashComparator());

        // Create cellStart array
        thrust::fill(d_cellStart.begin(), d_cellStart.end(), -1);
        const int threadsPerBlock = 256;
        const int numBlocks = (params.numParticles + threadsPerBlock - 1) / threadsPerBlock;
        find_cell_start_kernal<<<numBlocks, threadsPerBlock>>>(
            thrust::raw_pointer_cast(d_particleArray.data()),
            thrust::raw_pointer_cast(d_cellStart.data()),
            params.numParticles);
        
        // Parallel calculate denisty
        thrust::for_each(thrust::make_counting_iterator(0), thrust::make_counting_iterator(numParticles),
            CalcDensityFunctor(params, bounds, cellSize, gridRows, gridCols, d_cellStart, d_particleArray));

        // Parallel calculate pressure forces
        thrust::for_each(thrust::make_counting_iterator(0), thrust::make_counting_iterator(numParticles),
            CalcPressureFunctor(params, bounds, cellSize, gridRows, gridCols, d_cellStart, d_particleArray));

        // Apply other half of leapfrog
        thrust::for_each(d_particleArray.begin(), d_particleArray.end(),
            ApplyHalfGravityFunctor(params));

        // Bounce the particles
        thrust::for_each(d_particleArray.begin(), d_particleArray.end(),
            BounceParticleFunctor(params, bounds));

        // Stop execution time (dont count file saving)
        totalTime += std::chrono::duration_cast<mircoseconds>(timer::now() - start);

        // Save every saveFrames
        if (step % saveFrames == 0) {
            h_particleArray = d_particleArray;

            vector<float> positions(numParticles * numDimensions);
            vector<float> velocities(numParticles * numDimensions);
            for(int i = 0; i < numParticles; ++i) {
                positions[get_x_index(i)] = h_particleArray[i].x;
                positions[get_y_index(i)] = h_particleArray[i].y;
                velocities[get_x_index(i)] = h_particleArray[i].vx;
                velocities[get_y_index(i)] = h_particleArray[i].vy;
            }
            write_frame_to_gsd(handle, step, numParticles, bounds, positions, velocities, params);
        }
    }

    gsd_close(&handle);
    std::cout << "Total time(minus file saving) in 'us': " << totalTime.count() << std::endl;
    std::cout << "End of simulation" << std::endl;
    return 0;
}

// ==================== Custom CUDA Kernals ========================

// In parallel find the satrt of each cell in the list of sorted particles
__global__ void find_cell_start_kernal(const Particle* particleArray, int* cellStart, int numParticles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numParticles) return;
    int currentHash = particleArray[i].gridHash;
    if (i == 0) {
        cellStart[currentHash] = 0;
        return;
    }
    int previousHash = particleArray[i - 1].gridHash;
    if (currentHash != previousHash) {
        cellStart[currentHash] = i;
    }
}

// =================== Kernals ====================

// (r^2-d^2)^3 kernal
__device__ inline float ploy6_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float radius2 = radius * radius;
    float distance2 = dist * dist;
    float influence = radius2 - distance2;
    float norm = 4.0 / (M_PI * powf(radius, 8));
    return norm * (influence * influence * influence);
}

// (r-d)^3 kernal
__device__ inline float spikey_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float norm = 10.0 / (M_PI * powf(radius, 5));
    float influenece = radius - dist;
    return norm * (influenece * influenece * influenece);
}

// -3(r-d)^2
__device__ inline float spikey_kernal_derivative(float dist, float radius)
{
    if (dist >= radius || dist <= 0) return 0;
    float norm = 10.0 / (M_PI * powf(radius, 5));
    float influence = radius - dist;
    return -3.0 * norm * (influence * influence);
}

// Viscosity Kernal
__device__ inline float viscosity_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float norm = 45.0 / (M_PI * powf(radius, 6));
    float influence = radius - dist;
    return norm * influence;
}

// =============== helpers =======================

// For vectors
int get_x_index(int i)
{
    return i * numDimensions;
}

// For vectors
int get_y_index(int i)
{
    return i * numDimensions + 1;
}

__device__ inline float distance_from_two_points(float x1, float y1, float x2, float y2)
{
    float square_diff_x = (x2 - x1) * (x2 - x1);
    float square_diff_y = (y2 - y1) * (y2 - y1);
    float sum = square_diff_x + square_diff_y;
    float value = sqrtf(sum);
    return value;
}

// Simple calculation of pressure, depending on the target density
__device__ float density_to_pressure(float density, const SimulationParams &params)
{
    float densityError = density - params.targetDensity;
    float pressure = densityError * params.pressureMulti;
    return pressure;
}

// ============= Particle INITS ==================

// Randomly generates positions for each particle
void init_random_particle_pos_AoS(vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params)
{
    float radius = params.radius;
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> pos_x(bounds.xmin + radius, bounds.xmax - radius);
    std::uniform_real_distribution<float> pos_y(bounds.ymin + radius, bounds.ymax - radius);

    for (int i = 0; i < params.numParticles; i++) {
        particleArray[i].x = pos_x(gen);
        particleArray[i].y = pos_y(gen);

        particleArray[i].vx = 0.0;
        particleArray[i].vy = 0.0;

        particleArray[i].density = 0.0;
    }
}

void init_grid_particle_pos_AoS(thrust::host_vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params)
{
    int numParticles = params.numParticles;
    float radius = params.radius;
    float spacing = 2.0 * radius;

    float spawnWidth = (bounds.xmax - bounds.xmin) - (2 * radius);
    float spawnHeight = (bounds.ymax - bounds.ymin) - (2 * radius);

    if (spawnWidth < 0 || spawnHeight < 0) {
        std::cerr << "Bounding box is too small" << std::endl;
        exit(1);
    }

    int cols_fittable = static_cast<int>(floorf(spawnWidth / spacing)) + 1;
    int rows_fittable = static_cast<int>(floorf(spawnHeight / spacing)) + 1;
    int max_particles = cols_fittable * rows_fittable;

    if (numParticles > max_particles) {
        std::cerr << "Cannot fit " << numParticles << " particles in the box" << std::endl;
        std::cerr << "The maximum is " << max_particles << std::endl;
        exit(1);
    }

    int actual_cols = static_cast<int>(ceil(sqrt(numParticles)));
    int actual_rows = static_cast<int>(ceil(static_cast<float>(numParticles) / actual_cols));

    float block_width = (actual_cols - 1) * spacing;
    float block_height = (actual_rows - 1) * spacing;
    
    float offset_x = (spawnWidth - block_width) / 2.0f;
    float offset_y = (spawnHeight - block_height) / 2.0f;

    float startX = bounds.xmin + radius + offset_x;
    float startY = bounds.ymin + radius + offset_y;

    int currentParticle = 0;
    for (int i = 0; i < actual_rows && currentParticle < numParticles; ++i) {
        for (int j = 0; j < actual_cols && currentParticle < numParticles; ++j) {
            
            particleArray[currentParticle].x = startX + j * spacing;
            particleArray[currentParticle].y = startY + i * spacing;

            particleArray[currentParticle].vx = 0.0;
            particleArray[currentParticle].vy = 0.0;
            particleArray[currentParticle].density = 0.0;
            
            currentParticle++;
        }
    }
}

// ====================== Json ===================================

// Save metadata of simulation to JSON
void save_to_json(const std::string &filename, const SimulationParams &params, const box &bounds)
{
    cJSON *json = cJSON_CreateObject();
    cJSON_AddNumberToObject(json, "numParticles", params.numParticles);
    cJSON_AddNumberToObject(json, "timeStep", params.timeStep);
    cJSON_AddNumberToObject(json, "gravity", params.gravity);
    cJSON_AddNumberToObject(json, "elasticity", params.elasticity);
    cJSON_AddNumberToObject(json, "mass", params.mass);
    cJSON_AddNumberToObject(json, "radius", params.radius);
    cJSON_AddNumberToObject(json, "smoothingRadius", params.smoothingRadius);
    cJSON_AddNumberToObject(json, "targetDensity", params.targetDensity);
    cJSON_AddNumberToObject(json, "pressureMulti", params.pressureMulti);
    cJSON_AddNumberToObject(json, "xmax", bounds.xmax);
    cJSON_AddNumberToObject(json, "xmin", bounds.xmin);
    cJSON_AddNumberToObject(json, "ymax", bounds.ymax);
    cJSON_AddNumberToObject(json, "ymin", bounds.ymin);
    

    char *json_string = cJSON_Print(json);
    std::ofstream file(filename);
    if (file.is_open()) {
        file << json_string;
        file.close();
    } else {
        std::cerr << "Unable to save to JSON" << std::endl;
    }

    cJSON_Delete(json);
    free(json_string);
}

// ====================== GSD file formats =============================

// REF: Got chatGPT to make most of this, I removed the error handling (exessive)
void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<float> &positions, const vector<float> &velocities,
        const SimulationParams &params)
{
    int ret = 0;

    // The HOOMD GSD schema expects 3D data. We will pad the z-dimension with zeros.
    vector<float> pos_float(numParticles * 3, 0.0f);
    vector<float> vel_float(numParticles * 3, 0.0f);

    float box_center_x = (bounds.xmax + bounds.xmin) / 2.0f;
    float box_center_y = (bounds.ymax + bounds.ymin) / 2.0f;
    float box_center_z = (bounds.zmax + bounds.zmin) / 2.0f;

    for (uint32_t i = 0; i < numParticles; ++i) {
        // Store positions relative to the box center
        pos_float[i * 3 + 0] = static_cast<float>(positions[get_x_index(i)] - box_center_x);
        pos_float[i * 3 + 1] = static_cast<float>(positions[get_y_index(i)] - box_center_y);
        pos_float[i * 3 + 2] = 0.0f - box_center_z; // Assuming z is 0

        vel_float[i * 3 + 0] = static_cast<float>(velocities[get_x_index(i)]);
        vel_float[i * 3 + 1] = static_cast<float>(velocities[get_y_index(i)]);
    }

    // for saving the particle radius
    vector<float> diameters(numParticles, 2.0f * params.radius);

    // Config data
    float box_data[6] = {bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin, bounds.zmax - bounds.zmin, 0, 0, 0};
    ret = gsd_write_chunk(&handle, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, &box_data);
    ret = gsd_write_chunk(&handle, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);

    // Particle Data
    ret = gsd_write_chunk(&handle, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &numParticles);
    ret = gsd_write_chunk(&handle, "particles/position", GSD_TYPE_FLOAT, numParticles, 3, 0, pos_float.data());
    // ret = gsd_write_chunk(&handle, "particles/velocity", GSD_TYPE_FLOAT, numParticles, 3, 0, vel_float.data());
    ret = gsd_write_chunk(&handle, "particles/diameter", GSD_TYPE_FLOAT, numParticles, 1, 0, diameters.data());
    // Particle ID
    vector<uint32_t> typeids(numParticles, 0);
    ret = gsd_write_chunk(&handle, "particles/typeid", GSD_TYPE_UINT32, numParticles, 1, 0, typeids.data());
    
    // Particle types
    const char types[] = "A\0";
    ret = gsd_write_chunk(&handle, "particles/types", GSD_TYPE_INT8, 1, 2, 0, types);

    gsd_end_frame(&handle);
}
