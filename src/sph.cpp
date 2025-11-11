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


// Ensum to hold all version that can be run
enum Versions {
    SoA = 0, //n^2 approach
    SoA_LINKED = 1, // struct of arrays, with linkedlist grid
    AoS_LINKED = 2, // array of structs, with linkedlist grid
    AoS_SORTED = 3 // array of structs, with sorted grid
};

// ========================= Func Protos =======================================
void apply_gravity_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params);
void update_position_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params);
void bounce_particle_AoS(int i, vector<Particle> &particleArray, const box &bounds, const SimulationParams &params);

// hash + sorting
void build_grid_sorted(vector<Particle> &particleArray, vector<int> &cellStart, const SimulationParams &params,
    float cellSize, const box &bounds, int &gridRows, int &gridCols);
void calculate_density_sorted(int i, vector<Particle> &particleArray, const SimulationParams &params,
    vector<int> &cellStart, int &gridRows, int &gridCols, float cellSize, const box &bounds);
void calculate_pressure_sorted(int i, float &force_x, float &force_y, vector<Particle> &particleArray,
        const SimulationParams &params, const vector<int> &cellStart, int gridRows, int gridCols,
        float cellSize, const box &bounds);

// Kernals
float ploy6_kernal(float dist, float radius);
float spikey_kernal(float dist, float radius);
float spikey_kernal_derivative(float dist, float radius);
inline float viscosity_kernal(float dist, float radius);

// General
int get_x_index(int i);
int get_y_index(int i);
float distance_from_two_points(float x1, float y1, float x2, float y2);

// Inits
void init_random_particle_pos_AoS(vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params);
void init_grid_particle_pos_AoS(vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params);

// File saving
void save_to_json(const std::string &filename, const SimulationParams &params, const box &bounds);
void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<float> &positions, const vector<float> &velocities,
        const SimulationParams &params);

// ========================== Functions ========================================

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
        .numParticles = numParticles,
        .timeStep = timeStep,
        .gravity = gravity,
        .elasticity = elasticity,
        .mass = mass,
        .radius = radius,
        .smoothingRadius = smoothingRadius,
        .targetDensity = targetDensity,
        .pressureMulti = pressureMulti,
        .viscosity = viscosity
    };

    // 2 array linked list setup
    float cellSize = smoothingRadius;
    int gridRows = 0;
    int gridCols = 0;
    vector<int> head;
    vector<int> next;

    // setup for sorted
    vector<int> cellStart;

    // ========================== File Names ==============================
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
    // LinkedList grid + AoS approach
    vector<Particle> particleArray(numParticles);

    // init_random_particle_pos_AoS(particleArray, bounds, params); // random starting pos
    init_grid_particle_pos_AoS(particleArray, bounds, params);
    for (int step = 0; step < maxNumSteps; step++) {
        // Time execution
        auto start = timer::now();
    
        for (int i = 0; i < numParticles; i++) {
            // Apply gravity and look one step in the future
            particleArray[i].vy += params.gravity * params.timeStep * 0.5;
            update_position_AoS(i, particleArray, params);
        }

        // Create the linklist struct for grids/ cells
        build_grid_sorted(particleArray, cellStart, params, cellSize, bounds,
                gridRows, gridCols);

        for (int i = 0; i < numParticles; i++) {
            // Calc density
            calculate_density_sorted(i, particleArray, params, cellStart, gridRows,
                    gridCols, cellSize, bounds);
        }

        for (int i = 0; i < numParticles; i++) {
            // Calc pressure
            float forceX = 0.0;
            float forceY = 0.0;

            calculate_pressure_sorted(i, forceX, forceY, particleArray, params, cellStart,
                    gridRows, gridCols, cellSize, bounds);

            // a = f / density
            if (particleArray[i].density <= 0) continue;
            float accelX = forceX / particleArray[i].density;
            float accelY = forceY / particleArray[i].density;
            particleArray[i].vx += accelX * timeStep;
            particleArray[i].vy += accelY * timeStep;

        }

        for (int i = 0; i < numParticles; i++) {
            // Update postions and box collisions
            // update_position_AoS(i, particleArray, params);
            particleArray[i].vy += params.gravity * params.timeStep * 0.5;
            bounce_particle_AoS(i, particleArray, bounds, params);
        }

        // Stop execution time (dont count file saving)
        totalTime += std::chrono::duration_cast<mircoseconds>(timer::now() - start);

        // Save every saveFrames
        if (step % saveFrames == 0) {
            // write_particle_positions_to_csv_AoS(csvFile, particleArray, step);
            vector<float> positions(numParticles * numDimensions);
            vector<float> velocities(numParticles * numDimensions);
            for(int i = 0; i < numParticles; ++i) {
                positions[get_x_index(i)] = particleArray[i].x;
                positions[get_y_index(i)] = particleArray[i].y;
                velocities[get_x_index(i)] = particleArray[i].vx;
                velocities[get_y_index(i)] = particleArray[i].vy;
            }
            write_frame_to_gsd(handle, step, numParticles, bounds, positions, velocities, params);
        }
    }

    gsd_close(&handle);
    std::cout << "Total time(minus file saving) in 'us': " << totalTime.count() << std::endl;
    std::cout << "End of simulation" << std::endl;

    return 0;
}

// Simple calculation of pressure, depending on the target density
float density_to_pressure(float density, const SimulationParams &params)
{
    float densityError = density - params.targetDensity;
    float pressure = densityError * params.pressureMulti;
    return pressure;
}

// ============================= AoS implementation ============================

void apply_gravity_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params)
{
    particleArray[i].vy += params.gravity * params.timeStep;
}

void update_position_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params)
{
    particleArray[i].x += particleArray[i].vx * params.timeStep;
    particleArray[i].y += particleArray[i].vy * params.timeStep;
}

void bounce_particle_AoS(int i, vector<Particle> &particleArray, const box &bounds, const SimulationParams &params)
{
    // check walls in x
    if (particleArray[i].x - params.radius < bounds.xmin) {
        particleArray[i].x = bounds.xmin + params.radius;
        particleArray[i].vx *= -1 * params.elasticity;
    }
    if (particleArray[i].x + params.radius > bounds.xmax) {
        particleArray[i].x = bounds.xmax - params.radius;
        particleArray[i].vx *= -1 * params.elasticity;
    }
    // Check walls in y
    if (particleArray[i].y - params.radius < bounds.ymin) {
        particleArray[i].y = bounds.ymin + params.radius;
        particleArray[i].vy *= -1 * params.elasticity;
    }
    if (particleArray[i].y + params.radius > bounds.ymax) {
        particleArray[i].y = bounds.ymax - params.radius;
        particleArray[i].vy *= -1 * params.elasticity;
    }
}

// ======================== Hash and sort approach ======================================

void build_grid_sorted(vector<Particle> &particleArray, vector<int> &cellStart, const SimulationParams &params,
    float cellSize, const box &bounds, int &gridRows, int &gridCols)
{
    float width = bounds.xmax - bounds.xmin;
    float height = bounds.ymax - bounds.ymin;
    gridRows = std::max(1, (int)ceil(height / cellSize));
    gridCols = std::max(1, (int)ceil(width / cellSize));

    int numCells = gridRows * gridCols;
    cellStart.assign(numCells, -1);

    // Find hash for each particle
    for (int i = 0; i < params.numParticles; i++) {
        float x = particleArray[i].x;
        float y = particleArray[i].y;

        int cellX = (int)floor((x - bounds.xmin) / cellSize);
        int cellY = (int)floor((y - bounds.ymin) / cellSize);
        cellX = std::clamp(cellX, 0, gridCols - 1);
        cellY = std::clamp(cellY, 0, gridRows - 1);

        particleArray[i].gridHash = cellX + (cellY * gridCols);
    }

    // Sort now according to gridHash
    std::sort(particleArray.begin(), particleArray.end(), [](const Particle& a, const Particle& b) {
        return a.gridHash < b.gridHash;
    });

    // Find the start of each cell in the array of particles
    if (params.numParticles == 0) return;
    cellStart[particleArray[0].gridHash] = 0;
    for (int i = 1; i < params.numParticles; i++) {
        if (particleArray[i].gridHash != particleArray[i-1].gridHash) {
            cellStart[particleArray[i].gridHash] = i;
        }
    }
}

void calculate_density_sorted(int i, vector<Particle> &particleArray, const SimulationParams &params,
    vector<int> &cellStart, int &gridRows, int &gridCols, float cellSize, const box &bounds)
{
    float density = 0;

    float x1 = particleArray[i].x;
    float y1 = particleArray[i].y;

    // get current cell pos
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

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

void calculate_pressure_sorted(int i, float &force_x, float &force_y, vector<Particle> &particleArray,
        const SimulationParams &params, const vector<int> &cellStart, int gridRows, int gridCols,
        float cellSize, const box &bounds)
{
    force_x = 0.0;
    force_y = 0.0;

    float x1 = particleArray[i].x;
    float y1 = particleArray[i].y;

    // Get current particle velocity
    float vx1 = particleArray[i].vx;
    float vy1 = particleArray[i].vy;

    float density1 = particleArray[i].density;
    if (density1 <= 0.0) return;
    float pressure1 = density_to_pressure(density1, params);

    // get current cell pos
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

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
}

// =================== Kernals ====================

// (r^2-d^2)^3 kernal
inline float ploy6_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float radius2 = radius * radius;
    float distance2 = dist * dist;
    float influence = radius2 - distance2;
    float norm = 4.0 / (M_PI * pow(radius, 8));
    return norm * (influence * influence * influence);
}

// (r-d)^3 kernal
inline float spikey_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float norm = 10.0 / (M_PI * pow(radius, 5));
    float influenece = radius - dist;
    return norm * (influenece * influenece * influenece);
}

// -3(r-d)^2
inline float spikey_kernal_derivative(float dist, float radius)
{
    if (dist >= radius || dist <= 0) return 0;
    float norm = 10.0 / (M_PI * pow(radius, 5));
    float influence = radius - dist;
    return -3.0 * norm * (influence * influence);
}

// Viscosity Kernal
inline float viscosity_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float norm = 45.0 / (M_PI * pow(radius, 6));
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

inline float distance_from_two_points(float x1, float y1, float x2, float y2)
{
    float square_diff_x = (x2 - x1) * (x2 - x1);
    float square_diff_y = (y2 - y1) * (y2 - y1);
    float sum = square_diff_x + square_diff_y;
    float value = sqrt(sum);
    return value;
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

void init_grid_particle_pos_AoS(vector<Particle> &particleArray,
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

    int cols_fittable = static_cast<int>(floor(spawnWidth / spacing)) + 1;
    int rows_fittable = static_cast<int>(floor(spawnHeight / spacing)) + 1;
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
