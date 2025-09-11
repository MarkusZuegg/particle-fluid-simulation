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

// Struct to hold all varibales needed for a particle
// AoS implementation
typedef struct
{
    double x;
    double y;
    double vx;
    double vy;
    double density;
} Particle;


// ========================= Func Protos =======================================

void bounce_particle(int i, vector<double> &position, vector<double> &velocity,
        const box &bounds, float elasticity, float radius);
void init_particle_positions(vector<double> pos);
void update_particle_position(int i, vector<double> &position, vector<double> &velocity,
        float timeStep);
void apply_gravity(int i, vector<double> &velocity, float timeStep, float gravity);
void predicted_positions(int i, vector<double> &position, vector<double> &velocity,
        vector<double> &predictedPos, float timeStep);
void calculate_denisty(int i, vector<double> &densities, vector<double> &position,
        float mass, int numParticles, float radius);
float density_to_pressure(double density, float targetDensity, float pressureMulti);
void calculate_pressure(int i, double &force_x, double &force_y, vector<double> &postion,
                        vector<double> &densities, int numParticles, int mass, float radius,
                        float targetDensity, float pressureMulti);
void build_grid_linkedlist(vector<double> &position, vector<int> &head, vector<int> &next,
        int numParticles, float cellSize, const box &bounds, int &gridRows, int &gridCols);
void calculate_denisty_grid(int i, vector<double> &densities, vector<double> &position,
        float mass, float radius, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds);
void calculate_pressure_grid(int i, double &force_x, double &force_y, vector<double> &postion,
        vector<double> &densities, int mass, float radius,
        float targetDensity, float pressureMulti, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds);

int get_x_index(int i);
int get_y_index(int i);
float distance_from_two_points(double x1, double y1, double x2, double y2);

float ploy6_kernal(float dist, float radius);
float spikey_kernal(float dist, float radius);
float spikey_kernal_derivative(float dist, float radius);

void init_random_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float radius);

std::ofstream open_csv_file(const std::string& filename, int numParticles);
void write_particle_positions_to_csv(std::ofstream &file, const vector<double>& positions, const int step);
void save_to_json(const std::string &filename, const int numParticles, const float radius,
        const box &bounds);
void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<double> &positions, const vector<double> &velocities);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg

    // Total time
    mircoseconds totalTime{0};

    // ============= Simulation params ==================
    const int maxNumSteps = 1500;
    const int saveFrames = 1; // How many frames to save. != 0
    const int numParticles = 4096;
    const double timeStep = 1.0/120.0;

    const float gravity = 0;
    const float elasticity = 0.85;

    const int mass = 1; // Particle mass
    const float radius = 0.1; // for particles

    const float smoothingRadius = 1.5;
    const float targetDensity = 1.5;
    const float pressureMulti = 10;

    const box bounds = {40, 0, 40, 0, 0, 0}; //xmax, xmin, ymax, ymin, zmax, zmin

    float cellSize = smoothingRadius;
    int gridRows = 0;
    int gridCols = 0;
    vector<int> head;
    vector<int> next;

    // ===================== Vector arrays ================================
    vector<double> position(numParticles * numDimensions); // x,y,z interleaved * numParticles
    vector<double> velocity(numParticles * numDimensions); // vx,vy,vz interleaved * numParticles
    vector<double> densities(numParticles); // Densities storing for all particles
    vector<double> predictedPos(numParticles * numDimensions);

    // ========================== File Names ==============================
    std::string fileNameCsv = "logs/x" + std::string(argv[0]) + ".csv";
    std::string fileNameJson = "logs/x" + std::string(argv[0]) + ".json";
    std::string fileNameGsd = "logs/x" + std::string(argv[0]) + ".gsd";
    
    // ========================= CSV saving ===============================
    std::ofstream csvFile = open_csv_file(fileNameCsv, numParticles);
    save_to_json(fileNameJson, numParticles, radius, bounds);

    // Set the amount of digits saved from double
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
    init_random_particle_pos(position, velocity, numParticles, bounds, radius); // random starting pos

    // Setup which version of the code to run
    int original = 1;
    int linkedListGrid = 2;
    int linkedAoS = 3;

    for (int step = 0; step < maxNumSteps; step++) {
        // Time execution
        auto start = timer::now();
    
        for (int i = 0; i < numParticles; i++) {
            // Apply gravity and look one step in the future
            apply_gravity(i, velocity, timeStep, gravity);
            predicted_positions(i, position, velocity, predictedPos, timeStep);
        }

        // Create the linklist struct for grids/ cells
        build_grid_linkedlist(predictedPos, head, next, numParticles, cellSize, bounds,
                gridRows, gridCols);

        for (int i = 0; i < numParticles; i++) {
            // Calc density
            // calculate_denisty(i, densities, predictedPos, mass, numParticles, smoothingRadius);
            calculate_denisty_grid(i, densities, predictedPos, mass, smoothingRadius, head,
                    next, gridRows, gridCols, cellSize, bounds);
        }

        for (int i = 0; i < numParticles; i++) {
            // Calc pressure
            double forceX = 0.0;
            double forceY = 0.0;

            // calculate_pressure(i, forceX, forceY, predictedPos, densities,
            //         numParticles, mass, smoothingRadius, targetDensity, pressureMulti);
            calculate_pressure_grid(i, forceX, forceY, predictedPos, densities, mass,
                    smoothingRadius, targetDensity, pressureMulti, head, next, gridRows,
                    gridCols, cellSize, bounds);

            // a = f / density
            double accelX = forceX / densities[i];
            double accelY = forceY / densities[i];
            velocity[get_x_index(i)] += accelX * timeStep;
            velocity[get_y_index(i)] += accelY * timeStep;

        }

        for (int i = 0; i < numParticles; i++) {
            // Update postions and box collisions
            update_particle_position(i, position, velocity, timeStep);
            bounce_particle(i, position, velocity, bounds, elasticity, radius);
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

// ===================== SoA implementation ===================================

// Check if particle is out side box, if so bounce a particle off of the bounding box.
// Simple rebound by inverting the velocity and pushing the particle back inside. Can have elasticity if wanted.
// @param
// int i is the index of the particle to be checked/ bounced
void bounce_particle(int i, vector<double> &position, vector<double> &velocity,
        const box &bounds, float elasticity, float radius)
{
    // Get particle index for x and y
    int x = get_x_index(i);
    int y = get_y_index(i);

    // check walls in x
    if (position[x] - radius < bounds.xmin) {
        position[x] = bounds.xmin + radius;
        velocity[x] *= -1 * elasticity;
    }
    if (position[x] + radius > bounds.xmax) {
        position[x] = bounds.xmax - radius;
        velocity[x] *= -1 * elasticity;
    }
    // Check walls in y
    if (position[y] - radius < bounds.ymin) {
        position[y] = bounds.ymin + radius;
        velocity[y] *= -1 * elasticity;
    }
    if (position[y] + radius > bounds.ymax) {
        position[y] = bounds.ymax - radius;
        velocity[y] *= -1 * elasticity;
    }
}

void apply_gravity(int i, vector<double> &velocity, float timeStep, float gravity)
{
    int y = get_y_index(i);

    velocity[y] += gravity * timeStep; // vy += g * dt
}

// Updates the position of a particle depending on its velocity
void update_particle_position(int i, vector<double> &position, vector<double> &velocity,
        float timeStep)
{
    // Get particle index for x and y
    int x = get_x_index(i);
    int y = get_y_index(i);

    position[x] += velocity[x] * timeStep; // x += vx * dt
    position[y] += velocity[y] * timeStep; // y += vy * dt
}

// Just goes forward by one time step and returns the position the particle would be at.
// Stores this new positon in the predictedPos vector array
void predicted_positions(int i, vector<double> &position, vector<double> &velocity,
        vector<double> &predictedPos, float timeStep)
{
    int x = get_x_index(i);
    int y = get_y_index(i);

    predictedPos[x] = position[x] + (velocity[x] * timeStep);
    predictedPos[y] = position[y] + (velocity[y] * timeStep);
}

// Calculate the density for particle at index i.
// Goes through all particles in O(n^2), apply 
// @param i int index of particle
// @param densities vector array to store densities 
// @param position vector array for all particle positons
// @param mass int particle mass
// @param numParticle int number of particles
void calculate_denisty(int i, vector<double> &densities, vector<double> &position,
        float mass, int numParticles, float radius)
{
    double density = 0;

    double x1 = position[get_x_index(i)];
    double y1 = position[get_y_index(i)];


    // loop through all particles
    for (int j = 0; j < numParticles; j++) {
        double x2 = position[get_x_index(j)];
        double y2 = position[get_y_index(j)];

        float dist = distance_from_two_points(x1, y1, x2, y2);
        float influence = ploy6_kernal(dist, radius); // Kernal for density
        density += mass * influence;
    }
    densities[i] = density;
}

void calculate_pressure(int i, double &force_x, double &force_y, vector<double> &postion,
        vector<double> &densities, int numParticles, int mass, float radius,
        float targetDensity, float pressureMulti)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    double x1 = postion[get_x_index(i)];
    double y1 = postion[get_y_index(i)];

    // Get density and pressure for 1st particle
    double density1 = densities[i];
    float pressure1 = density_to_pressure(density1, targetDensity, pressureMulti);
    if (density1 <= 0.0) return;

    // Loop through all the other particles
    for (int j = 0; j < numParticles; j++) {
        if (i == j) continue;

        // get other particle postion
        double x2 = postion[get_x_index(j)];
        double y2 = postion[get_y_index(j)];

        // Get distance/ magnitiude
        float dist = distance_from_two_points(x1, y1, x2, y2);
        if (dist == 0) continue; // Either particle is same, or ontop of eachother

        // Get unit vectors
        double dirX = (x1 - x2) / dist;
        double dirY = (y1 - y2) / dist;

        // Get the slope of the kernal at the distance
        float slope = spikey_kernal_derivative(dist, radius); // Grab the gradient of pressure change
        // Get density and pressure for 2nd particle
        double density2 = densities[j];
        float pressure2 = density_to_pressure(density2, targetDensity, pressureMulti);

        float gradX = dirX * slope;
        float gradY = dirY * slope;
        // Symetric pressure calculation
        float coeff = -mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

        force_x += coeff * gradX;
        force_y += coeff * gradY;
    }
}

// Simple calculation of pressure, depending on the target density
float density_to_pressure(double density, float targetDensity, float pressureMulti)
{
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMulti;
    return pressure;
}

// ======================== Grid shenanigans ===================================
// REF: https://web.archive.org/web/20140725014123/https://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

void build_grid_linkedlist(vector<double> &position, vector<int> &head, vector<int> &next,
        int numParticles, float cellSize, const box &bounds, int &gridRows, int &gridCols)
{
    float width = bounds.xmax - bounds.xmin;
    float height = bounds.ymax - bounds.ymin;
    gridRows = std::max(1, (int)ceil(height / cellSize));
    gridCols = std::max(1, (int)ceil(width / cellSize));

    int numCells = gridRows * gridCols;
    // Reset head and next vector arrays
    head.assign(numCells, -1);
    next.assign(numParticles, -1);

    // Insert the particles
    for (int i = 0; i < numParticles; i++) {
        double x = position[get_x_index(i)];
        double y = position[get_y_index(i)];

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

void calculate_denisty_grid(int i, vector<double> &densities, vector<double> &position,
        float mass, float radius, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    double density = 0;

    double x1 = position[get_x_index(i)];
    double y1 = position[get_y_index(i)];

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
                double x2 = position[get_x_index(j)];
                double y2 = position[get_y_index(j)];

                float dist = distance_from_two_points(x1, y1, x2, y2);
                float influence = ploy6_kernal(dist, radius); // Kernal for density
                density += mass * influence;
            }
        }
    }
    densities[i] = density;
}

void calculate_pressure_grid(int i, double &force_x, double &force_y, vector<double> &postion,
        vector<double> &densities, int mass, float radius,
        float targetDensity, float pressureMulti, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    double x1 = postion[get_x_index(i)];
    double y1 = postion[get_y_index(i)];

    // Get density and pressure for 1st particle
    double density1 = densities[i];
    float pressure1 = density_to_pressure(density1, targetDensity, pressureMulti);
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
                double x2 = postion[get_x_index(j)];
                double y2 = postion[get_y_index(j)];

                // Get distance/ magnitiude
                float dist = distance_from_two_points(x1, y1, x2, y2);
                if (dist == 0) continue; // Either particle is same, or ontop of eachother

                // Get unit vectors
                double dirX = (x1 - x2) / dist;
                double dirY = (y1 - y2) / dist;

                // Get the slope of the kernal at the distance
                float slope = spikey_kernal_derivative(dist, radius); // Grab the gradient of pressure change
                // Get density and pressure for 2nd particle
                double density2 = densities[j];
                float pressure2 = density_to_pressure(density2, targetDensity, pressureMulti);

                float gradX = dirX * slope;
                float gradY = dirY * slope;
                // Symetric pressure calculation
                float coeff = -mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

                force_x += coeff * gradX;
                force_y += coeff * gradY;
            }
        }
    }
}

// ============================= AoS implementation ============================

void apply_gravity_AoS(Particle &particle, float timeStep, float gravity)
{
    particle.vy += gravity * timeStep;
}

void update_position_AoS(Particle &particle, float timeStep)
{
    particle.x = particle.vx * timeStep;
    particle.y = particle.vy * timeStep;
}

void bounce_particle_AoS(Particle &particle, const box &bounds, float elasticity, float radius)
{
    // check walls in x
    if (particle.x - radius < bounds.xmin) {
        particle.x = bounds.xmin + radius;
        particle.vx *= -1 * elasticity;
    }
    if (particle.x + radius > bounds.xmax) {
        particle.x = bounds.xmax - radius;
        particle.vx *= -1 * elasticity;
    }
    // Check walls in y
    if (particle.y - radius < bounds.ymin) {
        particle.y = bounds.ymin + radius;
        particle.vy *= -1 * elasticity;
    }
    if (particle.y + radius > bounds.ymax) {
        particle.y = bounds.ymax - radius;
        particle.vy *= -1 * elasticity;
    }
}

void build_grid_linkedlist_AoS(vector<Particle> &particleArray, vector<int> &head, vector<int> &next,
        int numParticles, float cellSize, const box &bounds, int &gridRows, int &gridCols)
{
    float width = bounds.xmax - bounds.xmin;
    float height = bounds.ymax - bounds.ymin;
    gridRows = std::max(1, (int)ceil(height / cellSize));
    gridCols = std::max(1, (int)ceil(width / cellSize));

    int numCells = gridRows * gridCols;
    // Reset head and next vector arrays
    head.assign(numCells, -1);
    next.assign(numParticles, -1);

    // Insert the particles
    for (int i = 0; i < numParticles; i++) {
        double x = particleArray[i].x;
        double y = particleArray[i].y;

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
        float mass, float radius, vector<int> &head, vector<int> &next,
        int gridRows, int gridCols, float cellSize, const box &bounds)
{
    double density = 0;

    double x1 = particleArray[i].x;
    double y1 = particleArray[i].y;

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
                double x2 = particleArray[j].x;
                double y2 = particleArray[j].y;

                float dist = distance_from_two_points(x1, y1, x2, y2);
                float influence = ploy6_kernal(dist, radius); // Kernal for density
                density += mass * influence;
            }
        }
    }
    particleArray[i].density = density;
}

void calculate_pressure_grid_Aos(int i, double &force_x, double &force_y, vector<Particle> &particleArray,
        int mass, float radius, float targetDensity, float pressureMulti,
        vector<int> &head, vector<int> &next, int gridRows, int gridCols,
        float cellSize, const box &bounds)
{
    force_x = 0.0;
    force_y = 0.0;

    // Get current particle position
    double x1 = particleArray[i].x;
    double y1 = particleArray[i].y;

    // Get density and pressure for 1st particle
    double density1 = particleArray[i].density;
    float pressure1 = density_to_pressure(density1, targetDensity, pressureMulti);
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
                double x2 = particleArray[j].x;
                double y2 = particleArray[j].y;

                // Get distance/ magnitiude
                float dist = distance_from_two_points(x1, y1, x2, y2);
                if (dist == 0) continue; // Either particle is same, or ontop of eachother

                // Get unit vectors
                double dirX = (x1 - x2) / dist;
                double dirY = (y1 - y2) / dist;

                // Get the slope of the kernal at the distance
                float slope = spikey_kernal_derivative(dist, radius); // Grab the gradient of pressure change
                // Get density and pressure for 2nd particle
                double density2 = particleArray[j].density;
                float pressure2 = density_to_pressure(density2, targetDensity, pressureMulti);

                float gradX = dirX * slope;
                float gradY = dirY * slope;
                // Symetric pressure calculation
                float coeff = -mass * (pressure1 / (density1 * density1) + pressure2 / (density2 * density2));

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

float distance_from_two_points(double x1, double y1, double x2, double y2)
{
    double square_diff_x = (x2 - x1) * (x2 - x1);
    double square_diff_y = (y2 - y1) * (y2 - y1);
    double sum = square_diff_x + square_diff_y;
    float value = sqrt(sum); //TODO make this quicker
    return value;
}

// Randomly generates positions for each particle
void init_random_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float radius)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> pos_x(bounds.xmin + radius, bounds.xmax - radius);
    std::uniform_real_distribution<double> pos_y(bounds.ymin + radius, bounds.ymax - radius);

    for (int i = 0; i < numParticles; i++) {
        int x = get_x_index(i);
        int y = get_y_index(i);
        
        position[x] = pos_x(gen);
        position[y] = pos_y(gen);

        velocity[x] = 0.0;
        velocity[y] = 0.0;
    }
}

// ========================== csv saving =======================================

std::ofstream open_csv_file(const std::string& filename, int numParticles)
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "step,";
        for (int i = 0; i < numParticles; i++) {
            file << "x" << i << ",y" << i << ","; // adds x1,y1,x2,y2 etc
        }
        file << '\n';
        return file;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return file;
    }
}

void write_particle_positions_to_csv(std::ofstream &file, const vector<double>& positions, const int step)
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

// ====================== Json ===================================

// Save metadata of simulation to JSON
void save_to_json(const std::string &filename, const int numParticles, const float radius,
        const box &bounds)
{
    cJSON *json = cJSON_CreateObject();
    cJSON_AddNumberToObject(json, "numParticles", numParticles);
    cJSON_AddNumberToObject(json, "radius", radius);
    cJSON_AddNumberToObject(json, "xmax", bounds.xmax);
    cJSON_AddNumberToObject(json, "xmin", bounds.xmin);
    cJSON_AddNumberToObject(json, "ymax", bounds.ymax);
    cJSON_AddNumberToObject(json, "ymin", bounds.ymin);
    cJSON_AddNumberToObject(json, "zmax", bounds.zmax);
    cJSON_AddNumberToObject(json, "zmin", bounds.zmin);

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
        const box &bounds, const vector<double> &positions, const vector<double> &velocities)
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
