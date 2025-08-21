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

using std::vector;

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
} box;

// ========================= Func Protos =======================================

void bounce_particle(int i, vector<double> &position, vector<double> &velocity,
        const box &bounds, float elasticity, float radius);
void init_particle_positions(vector<double> pos);
void update_particle_position(int i, vector<double> &position, vector<double> &velocity,
        float timeStep);
void init_random_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float radius);
void apply_gravity(int i, vector<double> &velocity, float timeStep, float gravity);

std::ofstream open_csv_file(const std::string& filename, int numParticles);
void write_particle_positions_to_csv(std::ofstream &file, const vector<double>& positions, const int step);
void save_to_json(const std::string &filename, const int numParticles, const float radius,
        const box &bounds);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg

    // ============= Simulation params ==================
    const int maxNumSteps = 500;
    const int saveFrames = 10; // How many frames to save. != 0
    const float gravity = -9.81;
    const float elasticity = 0.95;
    const int numParticles = 10;
    const float radius = 1;
    const double timeStep = 1.0/60.0;

    const box bounds = {100, 0, 100, 0}; //xmax, xmin, ymax, ymin

    // ===================== Vector arrays ================================
    vector<double> position(numParticles * 2); // x,y interleaved * numParticles
    vector<double> velocity(numParticles * 2); // vx,vy interleaved * numParticles
    
    // ========================= CSV saving ===============================
    std::string fileNameCsv = "logs/x" + std::string(argv[0]) + ".csv";
    std::string fileNameJson = "logs/x" + std::string(argv[0]) + ".json";
    std::ofstream csvFile = open_csv_file(fileNameCsv, numParticles);
    save_to_json(fileNameJson, numParticles, radius, bounds);

    // Set the amount of digits saved from double
    int sigfigs = 4;
    csvFile << std::fixed << std::setprecision(sigfigs);

    // ===================== Simulation loop ==================================
    init_random_particle_pos(position, velocity, numParticles, bounds, radius); // random starting pos
    for (int step = 0; step < maxNumSteps; step++) {
        for (int i = 0; i <numParticles; i++) {
            apply_gravity(i, velocity, timeStep, gravity);
            update_particle_position(i, position, velocity, timeStep);
            bounce_particle(i, position, velocity, bounds, elasticity, radius);
        }

        // Save every saveFrames
        if (step % saveFrames == 0) {
            write_particle_positions_to_csv(csvFile, position, step);
        }
    }
    
    csvFile.close();
    std::cout << "End of simulation" << std::endl;
    return 0;
}

// Check if particle is out side box, if so bounce a particle off of the bounding box.
// Simple rebound by inverting the velocity and pushing the particle back inside. Can have elasticity if wanted.
// @param
// int i is the index of the particle to be checked/ bounced
void bounce_particle(int i, vector<double> &position, vector<double> &velocity,
        const box &bounds, float elasticity, float radius)
{
    // Get particle index for x and y
    int x = 2 * i;
    int y = 2 * i + 1;

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
    // Get particle index for x and y
    // int x = 2 * i;
    int y = 2 * i + 1;

    velocity[y] += gravity * timeStep; // vy += g * dt
}

void update_particle_position(int i, vector<double> &position, vector<double> &velocity,
        float timeStep)
{
    // Get particle index for x and y
    int x = 2 * i;
    int y = 2 * i + 1;

    position[x] += velocity[x] * timeStep; // x += vx * dt
    position[y] += velocity[y] * timeStep; // y += vy * dt
}

void init_random_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float radius)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> pos_x(bounds.xmin + radius, bounds.xmax - radius);
    std::uniform_real_distribution<double> pos_y(bounds.ymin + radius, bounds.ymax - radius);

    for (int i = 0; i < numParticles; i++) {
        int x = 2*i;
        int y = 2*i + 1;
        
        position[x] = pos_x(gen);
        position[y] = pos_y(gen);

        velocity[x] = 0.0;
        velocity[y] = 0.0;
    }
}

// Smoothing kernel for calculating density
// Uses a simple x^3 kernel
float SmoothingKernel(float radius, float dist)
{
    float zero = 0.0;
    float value = std::max(zero, radius - dist);
    return value * value * value;
}

//=========================== csv saving =======================================

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

void save_to_json(const std::string &filename, const int numParticles, const float radius,
        const box &bounds)
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "{\n";
        file << "\"numParticles\":" << numParticles << ",\n";
        file << "\"radius\":" << radius << ",\n";
        file << "\"xmax\":" << bounds.xmax << ",\n";
        file << "\"xmin\":" << bounds.xmin << ",\n";
        file << "\"ymax\":" << bounds.ymax << ",\n";
        file << "\"ymin\":" << bounds.ymin << ",\n";
        file << "}\n";
        file.close();
    } else {
        std::cerr << "Cant save to json file" << std::endl;
    }
}
