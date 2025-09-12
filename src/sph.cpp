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
} Particle;

// Ensum to hold all version that can be run
enum Versions {
    ORGINAL = 0, //n^2 approach
    SoA_LINKED = 1, // struct of arrays, with linkedlist grid
    AoS_LINKED = 2 // array of structs, with linkedlist grid
};

// ========================= Func Protos =======================================

// SoA
void bounce_particle(int i, vector<float> &position, vector<float> &velocity,
        const box &bounds, const SimulationParams &params);
void apply_gravity(int i, vector<float> &velocity, const SimulationParams &params);
void update_particle_position(int i, vector<float> &position, vector<float> &velocity,
        const SimulationParams &params);
void predicted_positions(int i, vector<float> &position, vector<float> &velocity,
        vector<float> &predictedPos, const SimulationParams &params);
void calculate_denisty(int i, vector<float> &densities, vector<float> &position,
        const SimulationParams &params);
void calculate_pressure(int i, float &force_x, float &force_y, vector<float> &postion,
                        vector<float> &densities, const SimulationParams &params);
float density_to_pressure(float density, const SimulationParams &params);

// Linked list grid
void build_grid_linkedlist(vector<float> &position, vector<int> &head, vector<int> &next,
        const SimulationParams &params, float cellSize, const box &bounds, int &gridRows, int &gridCols);
void calculate_denisty_grid(int i, vector<float> &densities, vector<float> &position,
        const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds);
void calculate_pressure_grid(int i, float &force_x, float &force_y, vector<float> &postion,
        vector<float> &densities, const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds);

// AoS
void apply_gravity_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params);
void update_position_AoS(int i, vector<Particle> &particleArray, const SimulationParams &params);
void bounce_particle_AoS(int i, vector<Particle> &particleArray, const box &bounds, const SimulationParams &params);
void build_grid_linkedlist_AoS(vector<Particle> &particleArray, vector<int> &head, vector<int> &next,
        const SimulationParams &params, float cellSize, const box &bounds, int &gridRows, int &gridCols);
void calculate_denisty_grid_AoS(int i, vector<Particle> &particleArray,
        const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds);
void calculate_pressure_grid_AoS(int i, float &force_x, float &force_y, vector<Particle> &particleArray,
        const SimulationParams &params, vector<int> &head, vector<int> &next, int gridRows,
        int gridCols, float cellSize, const box &bounds);

// Kernals
float ploy6_kernal(float dist, float radius);
float spikey_kernal(float dist, float radius);
float spikey_kernal_derivative(float dist, float radius);

// General
int get_x_index(int i);
int get_y_index(int i);
float distance_from_two_points(float x1, float y1, float x2, float y2);
inline float fast_sqrt(float x);
float fast_distance(float x1, float y1, float x2, float y2);

// Inits
void init_random_particle_pos(vector<float> &position, vector<float> &velocity,
        const box &bounds, const SimulationParams &params);
void init_random_particle_pos_AoS(vector<Particle> &particleArray,
        const box &bounds, const SimulationParams &params);

// File saving
std::ofstream open_csv_file(const std::string& filename, const SimulationParams &params);
void write_particle_positions_to_csv(std::ofstream &file, const vector<float>& positions, const int step);
void write_particle_positions_to_csv_AoS(std::ofstream &file, const vector<Particle> &particleArray, const int step);
void save_to_json(const std::string &filename, const SimulationParams &params, const box &bounds);
void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<float> &positions, const vector<float> &velocities);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg

    // Total time
    mircoseconds totalTime{0};

    // ============= Simulation params ==================
    const int maxNumSteps = 240;
    const int saveFrames = 1; // How many frames to save. != 0
    const int numParticles = 8000;
    const float timeStep = 1.0/120.0;

    const float gravity = 0;
    const float elasticity = 0.85;

    const int mass = 1; // Particle mass
    const float radius = 0.1; // for particles

    const float smoothingRadius = 1.5;
    const float targetDensity = 1.5;
    const float pressureMulti = 10;

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
        .pressureMulti = pressureMulti
    };

    const box bounds = {40, 0, 40, 0, 0, 0}; //xmax, xmin, ymax, ymin, zmax, zmin

    // 2 array linked list setup
    float cellSize = smoothingRadius;
    int gridRows = 0;
    int gridCols = 0;
    vector<int> head;
    vector<int> next;

    // ========================== File Names ==============================
    std::string fileNameCsv = "logs/x" + std::string(argv[0]) + ".csv";
    std::string fileNameJson = "logs/x" + std::string(argv[0]) + ".json";
    std::string fileNameGsd = "logs/x" + std::string(argv[0]) + ".gsd";
    
    // ========================= CSV saving ===============================
    std::ofstream csvFile = open_csv_file(fileNameCsv, params);
    save_to_json(fileNameJson, params, bounds);

    // Set the amount of digits saved from float
    int sigfigs = 4;
    csvFile << std::fixed << std::setprecision(sigfigs);

    // =========================== GSD ====================================
    gsd_handle handle;
    int ret = gsd_create_and_open(&handle, fileNameGsd.c_str(), "SPH simulation",
            "hoomd", gsd_make_version(1, 0), GSD_OPEN_READWRITE, 0);
    if (ret != 0) {
        std::cerr << "Error opening GSD file" << std::endl;
        return -1;
    }

    // ===================== Simulation loop ==================================

    // Change the line below to switch between each version
    int runType = AoS_LINKED;

    // LinkedList grid + SoA approach
    if (runType == SoA_LINKED) {
        // ===================== Vector arrays ================================
        vector<float> position(numParticles * numDimensions); // x,y,z interleaved * numParticles
        vector<float> velocity(numParticles * numDimensions); // vx,vy,vz interleaved * numParticles
        vector<float> densities(numParticles); // Densities storing for all particles
        vector<float> predictedPos(numParticles * numDimensions);

        init_random_particle_pos(position, velocity, bounds, params); // random starting pos
        for (int step = 0; step < maxNumSteps; step++) {
            // Time execution
            auto start = timer::now();
        
            for (int i = 0; i < numParticles; i++) {
                // Apply gravity and look one step in the future
                apply_gravity(i, velocity, params);
                predicted_positions(i, position, velocity, predictedPos, params);
            }

            // Create the linklist struct for grids/ cells
            build_grid_linkedlist(predictedPos, head, next, params, cellSize, bounds,
                    gridRows, gridCols);

            for (int i = 0; i < numParticles; i++) {
                // Calc density
                // calculate_denisty(i, densities, predictedPos, mass, numParticles, smoothingRadius);
                calculate_denisty_grid(i, densities, predictedPos, params, head,
                        next, gridRows, gridCols, cellSize, bounds);
            }

            for (int i = 0; i < numParticles; i++) {
                // Calc pressure
                float forceX = 0.0;
                float forceY = 0.0;

                // calculate_pressure(i, forceX, forceY, predictedPos, densities,
                //         numParticles, mass, smoothingRadius, targetDensity, pressureMulti);
                calculate_pressure_grid(i, forceX, forceY, predictedPos, densities, params,
                        head, next, gridRows, gridCols, cellSize, bounds);

                // a = f / density
                float accelX = forceX / densities[i];
                float accelY = forceY / densities[i];
                velocity[get_x_index(i)] += accelX * timeStep;
                velocity[get_y_index(i)] += accelY * timeStep;

            }

            for (int i = 0; i < numParticles; i++) {
                // Update postions and box collisions
                update_particle_position(i, position, velocity, params);
                bounce_particle(i, position, velocity, bounds, params);
            }

            // Stop execution time (dont count file saving)
            totalTime += std::chrono::duration_cast<mircoseconds>(timer::now() - start);

            // Save every saveFrames
            if (step % saveFrames == 0) {
                write_particle_positions_to_csv(csvFile, position, step);
                write_frame_to_gsd(handle, step, numParticles, bounds, position, velocity);
            }
        }
        
        csvFile.close();
        gsd_close(&handle);
        std::cout << "Total time(minus file saving) in 'us': " << totalTime.count() << std::endl;
        std::cout << "End of simulation" << std::endl;
        return 0;
    }

    // LinkedList grid + AoS approach
    if (runType == AoS_LINKED) {
        vector<Particle> particleArray(numParticles);

        init_random_particle_pos_AoS(particleArray, bounds, params); // random starting pos
        for (int step = 0; step < maxNumSteps; step++) {
            // Time execution
            auto start = timer::now();
        
            for (int i = 0; i < numParticles; i++) {
                // Apply gravity and look one step in the future
                apply_gravity_AoS(i, particleArray, params);
                // predicted_positions(i, particleArray, predictedPos, params);
            }

            // Create the linklist struct for grids/ cells
            build_grid_linkedlist_AoS(particleArray, head, next, params, cellSize, bounds,
                    gridRows, gridCols);

            for (int i = 0; i < numParticles; i++) {
                // Calc density
                // calculate_denisty(i, densities, predictedPos, mass, numParticles, smoothingRadius);
                calculate_denisty_grid_AoS(i, particleArray, params, head,
                        next, gridRows, gridCols, cellSize, bounds);
            }

            for (int i = 0; i < numParticles; i++) {
                // Calc pressure
                float forceX = 0.0;
                float forceY = 0.0;

                // calculate_pressure(i, forceX, forceY, predictedPos, densities,
                //         numParticles, mass, smoothingRadius, targetDensity, pressureMulti);
                calculate_pressure_grid_AoS(i, forceX, forceY, particleArray, params,
                        head, next, gridRows, gridCols, cellSize, bounds);

                // a = f / density
                float accelX = forceX / particleArray[i].density;
                float accelY = forceY / particleArray[i].density;
                particleArray[i].vx += accelX * timeStep;
                particleArray[i].vy += accelY * timeStep;

            }

            for (int i = 0; i < numParticles; i++) {
                // Update postions and box collisions
                update_position_AoS(i, particleArray, params);
                bounce_particle_AoS(i, particleArray, bounds, params);
            }

            // Stop execution time (dont count file saving)
            totalTime += std::chrono::duration_cast<mircoseconds>(timer::now() - start);

            // Save every saveFrames
            if (step % saveFrames == 0) {
                write_particle_positions_to_csv_AoS(csvFile, particleArray, step);
                // write_frame_to_gsd(handle, step, numParticles, bounds, position, velocity);
            }
        }
        
        csvFile.close();
        gsd_close(&handle);
        std::cout << "Total time(minus file saving) in 'us': " << totalTime.count() << std::endl;
        std::cout << "End of simulation" << std::endl;
        return 0;

    }

    // N^2 approach + SoA
    if (runType == ORGINAL) {
        // ===================== Vector arrays ================================
        vector<float> position(numParticles * numDimensions); // x,y,z interleaved * numParticles
        vector<float> velocity(numParticles * numDimensions); // vx,vy,vz interleaved * numParticles
        vector<float> densities(numParticles); // Densities storing for all particles
        vector<float> predictedPos(numParticles * numDimensions);

        init_random_particle_pos(position, velocity, bounds, params); // random starting pos
        for (int step = 0; step < maxNumSteps; step++) {
            // Time execution
            auto start = timer::now();
        
            for (int i = 0; i < numParticles; i++) {
                // Apply gravity and look one step in the future
                apply_gravity(i, velocity, params);
                predicted_positions(i, position, velocity, predictedPos, params);
            }

            for (int i = 0; i < numParticles; i++) {
                // Calc density
                calculate_denisty(i, densities, predictedPos, params);
            }

            for (int i = 0; i < numParticles; i++) {
                // Calc pressure
                float forceX = 0.0;
                float forceY = 0.0;

                calculate_pressure(i, forceX, forceY, predictedPos, densities, params);

                // a = f / density
                float accelX = forceX / densities[i];
                float accelY = forceY / densities[i];
                velocity[get_x_index(i)] += accelX * timeStep;
                velocity[get_y_index(i)] += accelY * timeStep;

            }

            for (int i = 0; i < numParticles; i++) {
                // Update postions and box collisions
                update_particle_position(i, position, velocity, params);
                bounce_particle(i, position, velocity, bounds, params);
            }

            // Stop execution time (dont count file saving)
            totalTime += std::chrono::duration_cast<mircoseconds>(timer::now() - start);

            // Save every saveFrames
            if (step % saveFrames == 0) {
                write_particle_positions_to_csv(csvFile, position, step);
                write_frame_to_gsd(handle, step, numParticles, bounds, position, velocity);
            }
        }
        
        csvFile.close();
        gsd_close(&handle);
        std::cout << "Total time(minus file saving) in 'us': " << totalTime.count() << std::endl;
        std::cout << "End of simulation" << std::endl;
        return 0;
    }

    // Didnt go into the if statements/ versions
    std::cout << "Didnt pick a valid version" << std::endl;
    return 0;
}

// ===================== SoA implementation ===================================

// Check if particle is out side box, if so bounce a particle off of the bounding box.
// Simple rebound by inverting the velocity and pushing the particle back inside. Can have elasticity if wanted.
// @param
// int i is the index of the particle to be checked/ bounced
void bounce_particle(int i, vector<float> &position, vector<float> &velocity,
        const box &bounds, const SimulationParams &params)
{
    // Get particle index for x and y
    int x = get_x_index(i);
    int y = get_y_index(i);

    // check walls in x
    if (position[x] - params.radius < bounds.xmin) {
        position[x] = bounds.xmin + params.radius;
        velocity[x] *= -1 * params.elasticity;
    }
    if (position[x] + params.radius > bounds.xmax) {
        position[x] = bounds.xmax - params.radius;
        velocity[x] *= -1 * params.elasticity;
    }
    // Check walls in y
    if (position[y] - params.radius < bounds.ymin) {
        position[y] = bounds.ymin + params.radius;
        velocity[y] *= -1 * params.elasticity;
    }
    if (position[y] + params.radius > bounds.ymax) {
        position[y] = bounds.ymax - params.radius;
        velocity[y] *= -1 * params.elasticity;
    }
}

void apply_gravity(int i, vector<float> &velocity, const SimulationParams &params)
{
    int y = get_y_index(i);

    velocity[y] += params.gravity * params.timeStep; // vy += g * dt
}

// Updates the position of a particle depending on its velocity
void update_particle_position(int i, vector<float> &position, vector<float> &velocity,
        const SimulationParams &params)
{
    // Get particle index for x and y
    int x = get_x_index(i);
    int y = get_y_index(i);

    position[x] += velocity[x] * params.timeStep; // x += vx * dt
    position[y] += velocity[y] * params.timeStep; // y += vy * dt
}

// Just goes forward by one time step and returns the position the particle would be at.
// Stores this new positon in the predictedPos vector array
void predicted_positions(int i, vector<float> &position, vector<float> &velocity,
        vector<float> &predictedPos, const SimulationParams &params)
{
    int x = get_x_index(i);
    int y = get_y_index(i);

    predictedPos[x] = position[x] + (velocity[x] * params.timeStep);
    predictedPos[y] = position[y] + (velocity[y] * params.timeStep);
}

// Calculate the density for particle at index i.
// Goes through all particles in O(n^2), apply 
// @param i int index of particle
// @param densities vector array to store densities 
// @param position vector array for all particle positons
// @param mass int particle mass
// @param numParticle int number of particles
void calculate_denisty(int i, vector<float> &densities, vector<float> &position,
        const SimulationParams &params)
{
    float density = 0;

    float x1 = position[get_x_index(i)];
    float y1 = position[get_y_index(i)];


    // loop through all particles
    for (int j = 0; j < params.numParticles; j++) {
        float x2 = position[get_x_index(j)];
        float y2 = position[get_y_index(j)];

        float dist = distance_from_two_points(x1, y1, x2, y2);
        float influence = ploy6_kernal(dist, params.smoothingRadius); // Kernal for density
        density += params.mass * influence;
    }
    densities[i] = density;
}

void calculate_pressure(int i, float &force_x, float &force_y, vector<float> &postion,
                        vector<float> &densities, const SimulationParams &params)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    float x1 = postion[get_x_index(i)];
    float y1 = postion[get_y_index(i)];

    // Get density and pressure for 1st particle
    float density1 = densities[i];
    float pressure1 = density_to_pressure(density1, params);
    if (density1 <= 0.0) return;

    // Loop through all the other particles
    for (int j = 0; j < params.numParticles; j++) {
        if (i == j) continue;

        // get other particle postion
        float x2 = postion[get_x_index(j)];
        float y2 = postion[get_y_index(j)];

        // Get distance/ magnitiude
        float dist = distance_from_two_points(x1, y1, x2, y2);
        if (dist == 0) continue; // Either particle is same, or ontop of eachother

        // Get unit vectors
        float dirX = (x1 - x2) / dist;
        float dirY = (y1 - y2) / dist;

        // Get the slope of the kernal at the distance
        float slope = spikey_kernal_derivative(dist, params.smoothingRadius); // Grab the gradient of pressure change
        // Get density and pressure for 2nd particle
        float density2 = densities[j];
        float pressure2 = density_to_pressure(density2, params);

        float gradX = dirX * slope;
        float gradY = dirY * slope;
        // Symetric pressure calculation
        float coeff = -params.mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

        force_x += coeff * gradX;
        force_y += coeff * gradY;
    }
}

// Simple calculation of pressure, depending on the target density
float density_to_pressure(float density, const SimulationParams &params)
{
    float densityError = density - params.targetDensity;
    float pressure = densityError * params.pressureMulti;
    return pressure;
}

// ======================== Grid shenanigans ===================================
// REF: https://web.archive.org/web/20140725014123/https://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

void build_grid_linkedlist(vector<float> &position, vector<int> &head, vector<int> &next,
        const SimulationParams &params, float cellSize, const box &bounds, int &gridRows, int &gridCols)
{
    float width = bounds.xmax - bounds.xmin;
    float height = bounds.ymax - bounds.ymin;
    gridRows = std::max(1, (int)ceil(height / cellSize));
    gridCols = std::max(1, (int)ceil(width / cellSize));

    int numCells = gridRows * gridCols;
    // Reset head and next vector arrays
    head.assign(numCells, -1);
    next.assign(params.numParticles, -1);

    // Insert the particles
    for (int i = 0; i < params.numParticles; i++) {
        float x = position[get_x_index(i)];
        float y = position[get_y_index(i)];

        int cellX = (int)floor((x - bounds.xmin) / cellSize);
        int cellY = (int)floor((y - bounds.ymin) / cellSize);
        cellX = std::clamp(cellX, 0, gridCols - 1);
        cellY = std::clamp(cellY, 0, gridRows - 1);

        int cellIndex = cellX + (cellY * gridCols);
        // Swap head with new particle (do linkedlist things)
        next[i] = head[cellIndex];
        head[cellIndex] = i;
    }
}

void calculate_denisty_grid(int i, vector<float> &densities, vector<float> &position,
        const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    float density = 0;

    float x1 = position[get_x_index(i)];
    float y1 = position[get_y_index(i)];

    // Get first particle cell position
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

    // Go through a 3x3 of cells around the particle cell
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int newCellX = cellX + dx;
            int newCellY = cellY + dy;
            if (newCellX < 0 || newCellX >= gridCols) continue; // check x is inbounds
            if (newCellY < 0 || newCellY >= gridRows) continue; // check y

            int neighbourCell = newCellX + (newCellY * gridCols);
            // Now loop through all particles in this cell
            for (int j = head[neighbourCell]; j != -1; j = next[j]) {
                float x2 = position[get_x_index(j)];
                float y2 = position[get_y_index(j)];

                float dist = distance_from_two_points(x1, y1, x2, y2);
                float influence = ploy6_kernal(dist, params.smoothingRadius); // Kernal for density
                density += params.mass * influence;
            }
        }
    }
    densities[i] = density;
}

void calculate_pressure_grid(int i, float &force_x, float &force_y, vector<float> &postion,
        vector<float> &densities, const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    float x1 = postion[get_x_index(i)];
    float y1 = postion[get_y_index(i)];

    // Get density and pressure for 1st particle
    float density1 = densities[i];
    float pressure1 = density_to_pressure(density1, params);
    if (density1 <= 0.0) return;

    // Get first particle cell position
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

    // Go through a 3x3 of cells around the particle cell
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int newCellX = cellX + dx;
            int newCellY = cellY + dy;
            if (newCellX < 0 || newCellX >= gridCols) continue; // check x is inbounds
            if (newCellY < 0 || newCellY >= gridRows) continue; // check y

            int neighbourCell = newCellX + (newCellY * gridCols);
            // Now loop through all particles in this cell
            for (int j = head[neighbourCell]; j != -1; j = next[j]) {
                if (j == i) continue;

                // get other particle postion
                float x2 = postion[get_x_index(j)];
                float y2 = postion[get_y_index(j)];

                // Get distance/ magnitiude
                float dist = distance_from_two_points(x1, y1, x2, y2);
                if (dist == 0) continue; // Either particle is same, or ontop of eachother

                // Get unit vectors
                float dirX = (x1 - x2) / dist;
                float dirY = (y1 - y2) / dist;

                // Get the slope of the kernal at the distance
                float slope = spikey_kernal_derivative(dist, params.smoothingRadius); // Grab the gradient of pressure change
                // Get density and pressure for 2nd particle
                float density2 = densities[j];
                float pressure2 = density_to_pressure(density2, params);

                float gradX = dirX * slope;
                float gradY = dirY * slope;
                // Symetric pressure calculation
                float coeff = -params.mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

                force_x += coeff * gradX;
                force_y += coeff * gradY;
            }
        }
    }
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

void build_grid_linkedlist_AoS(vector<Particle> &particleArray, vector<int> &head, vector<int> &next,
        const SimulationParams &params, float cellSize, const box &bounds, int &gridRows, int &gridCols)
{
    float width = bounds.xmax - bounds.xmin;
    float height = bounds.ymax - bounds.ymin;
    gridRows = std::max(1, (int)ceil(height / cellSize));
    gridCols = std::max(1, (int)ceil(width / cellSize));

    int numCells = gridRows * gridCols;
    // Reset head and next vector arrays
    head.assign(numCells, -1);
    next.assign(params.numParticles, -1);

    // Insert the particles
    for (int i = 0; i < params.numParticles; i++) {
        float x = particleArray[i].x;
        float y = particleArray[i].y;

        int cellX = (int)floor((x - bounds.xmin) / cellSize);
        int cellY = (int)floor((y - bounds.ymin) / cellSize);
        cellX = std::clamp(cellX, 0, gridCols - 1);
        cellY = std::clamp(cellY, 0, gridRows - 1);

        int cellIndex = cellX + (cellY * gridCols);
        // Swap head with new particle (do linkedlist things)
        next[i] = head[cellIndex];
        head[cellIndex] = i;
    }
}

void calculate_denisty_grid_AoS(int i, vector<Particle> &particleArray,
        const SimulationParams &params, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    float density = 0;

    float x1 = particleArray[i].x;
    float y1 = particleArray[i].y;

    // Get first particle cell position
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

    // Go through a 3x3 of cells around the particle cell
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int newCellX = cellX + dx;
            int newCellY = cellY + dy;
            if (newCellX < 0 || newCellX >= gridCols) continue; // check x is inbounds
            if (newCellY < 0 || newCellY >= gridRows) continue; // check y

            int neighbourCell = newCellX + (newCellY * gridCols);
            // Now loop through all particles in this cell
            for (int j = head[neighbourCell]; j != -1; j = next[j]) {
                float x2 = particleArray[j].x;
                float y2 = particleArray[j].y;

                float dist = fast_distance(x1, y1, x2, y2);
                float influence = ploy6_kernal(dist, params.smoothingRadius); // Kernal for density
                density += params.mass * influence;
            }
        }
    }
    particleArray[i].density = density;
}

void calculate_pressure_grid_AoS(int i, float &force_x, float &force_y, vector<Particle> &particleArray,
        const SimulationParams &params, vector<int> &head, vector<int> &next, int gridRows,
        int gridCols, float cellSize, const box &bounds)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    float x1 = particleArray[i].x;
    float y1 = particleArray[i].y;

    // Get density and pressure for 1st particle
    float density1 = particleArray[i].density;
    float pressure1 = density_to_pressure(density1, params);
    if (density1 <= 0.0) return;

    // Get first particle cell position
    int cellX = (int)floor((x1 - bounds.xmin) / cellSize);
    int cellY = (int)floor((y1 - bounds.ymin) / cellSize);
    cellX = std::clamp(cellX, 0, gridCols - 1);
    cellY = std::clamp(cellY, 0, gridRows - 1);

    // Go through a 3x3 of cells around the particle cell
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int newCellX = cellX + dx;
            int newCellY = cellY + dy;
            if (newCellX < 0 || newCellX >= gridCols) continue; // check x is inbounds
            if (newCellY < 0 || newCellY >= gridRows) continue; // check y

            int neighbourCell = newCellX + (newCellY * gridCols);
            // Now loop through all particles in this cell
            for (int j = head[neighbourCell]; j != -1; j = next[j]) {
                if (j == i) continue;

                // get other particle postion
                float x2 = particleArray[j].x;
                float y2 = particleArray[j].y;

                // Get distance/ magnitiude
                float dist = fast_distance(x1, y1, x2, y2);
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
            }
        }
    }
}

// =================== Kernals ====================

// (r^2-d^2)^3 kernal
float ploy6_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float radius2 = radius * radius;
    float distance2 = dist * dist;
    float influence = radius2 - distance2;
    float norm = 4.0 / (M_PI * pow(radius, 8));
    return norm * (influence * influence * influence);
}

// (r-d)^3 kernal
float spikey_kernal(float dist, float radius)
{
    if (dist >= radius) return 0;
    float norm = 10.0 / (M_PI * pow(radius, 5));
    float influenece = radius - dist;
    return norm * (influenece * influenece * influenece);
}

// -3(r-d)^2
float spikey_kernal_derivative(float dist, float radius)
{
    if (dist >= radius || dist <= 0) return 0;
    float norm = 10.0 / (M_PI * pow(radius, 5));
    float influence = radius - dist;
    return -3.0 * norm * (influence * influence);
}

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

float distance_from_two_points(float x1, float y1, float x2, float y2)
{
    float square_diff_x = (x2 - x1) * (x2 - x1);
    float square_diff_y = (y2 - y1) * (y2 - y1);
    float sum = square_diff_x + square_diff_y;
    float value = sqrt(sum); //TODO make this quicker
    return value;
}

// Fast invverse sqrt from Quake 3 baby!!!!
inline float fast_sqrt(float x)
{
    float xhalf = 0.5f * x;
    int i = *(int*)&x;               // evil bit-level hack
    i = 0x5f3759df - (i >> 1);       // magic constant
    x = *(float*)&i;
    x = x * (1.5f - xhalf * x * x);  // 1st Newton iteration
    return 1.0f / x;
}

float fast_distance(float x1, float y1, float x2, float y2)
{
    float dx = x2 - x1;
    float dy = y2 - y1;
    float sum = dx*dx + dy*dy;
    return fast_sqrt(sum);
}

// Randomly generates positions for each particle
void init_random_particle_pos(vector<float> &position, vector<float> &velocity,
        const box &bounds, const SimulationParams &params)
{
    float radius = params.radius;
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> pos_x(bounds.xmin + radius, bounds.xmax - radius);
    std::uniform_real_distribution<float> pos_y(bounds.ymin + radius, bounds.ymax - radius);

    for (int i = 0; i < params.numParticles; i++) {
        int x = get_x_index(i);
        int y = get_y_index(i);
        
        position[x] = pos_x(gen);
        position[y] = pos_y(gen);

        velocity[x] = 0.0;
        velocity[y] = 0.0;
    }
}

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

// ========================== csv saving =======================================

std::ofstream open_csv_file(const std::string& filename, const SimulationParams &params)
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "step,";
        for (int i = 0; i < params.numParticles; i++) {
            file << "x" << i << ",y" << i << ","; // adds x1,y1,x2,y2 etc
        }
        file << '\n';
        return file;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return file;
    }
}

void write_particle_positions_to_csv(std::ofstream &file, const vector<float>& positions, const int step)
{
    if (file.is_open()) {
        file << step << ',';
        for (size_t i = 0; i < positions.size(); ++i) {
            file << positions[i] << (i == positions.size() - 1 ? "" : ",");
        }
        file << "\n";
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void write_particle_positions_to_csv_AoS(std::ofstream &file, const vector<Particle> &particleArray, const int step)
{
    if (file.is_open()) {
        file << step << ',';
        for (size_t i = 0; i < particleArray.size(); ++i) {
            file << particleArray[i].x << particleArray[i].y << (i == particleArray.size() - 1 ? "" : ",");
        }
        file << "\n";
    } else {
        std::cerr << "Unable to open file" << std::endl;
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
        const box &bounds, const vector<float> &positions, const vector<float> &velocities)
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

    // Config data
    float box_data[6] = {bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin, bounds.zmax - bounds.zmin, 0, 0, 0};
    ret = gsd_write_chunk(&handle, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, &box_data);
    ret = gsd_write_chunk(&handle, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);

    // Particle Data
    ret = gsd_write_chunk(&handle, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &numParticles);
    ret = gsd_write_chunk(&handle, "particles/position", GSD_TYPE_FLOAT, numParticles, 3, 0, pos_float.data());
    ret = gsd_write_chunk(&handle, "particles/velocity", GSD_TYPE_FLOAT, numParticles, 3, 0, vel_float.data());

    // Particle ID
    vector<uint32_t> typeids(numParticles, 0);
    ret = gsd_write_chunk(&handle, "particles/typeid", GSD_TYPE_UINT32, numParticles, 1, 0, typeids.data());
    
    // Particle types
    const char types[] = "A\0";
    ret = gsd_write_chunk(&handle, "particles/types", GSD_TYPE_INT8, 1, 2, 0, types);

    gsd_end_frame(&handle);
}
