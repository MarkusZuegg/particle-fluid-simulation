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

void write_frame_to_gsd(gsd_handle &handle, uint32_t step, uint32_t numParticles,
        const box &bounds, const vector<double> &positions, const vector<double> &velocities);

int get_x_index(int i);
int get_y_index(int i);
int get_z_index(int i);

float smoothing_kernel(float dist, float radius);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg

    // ============= Simulation params ==================
    const int maxNumSteps = 500;
    const int saveFrames = 1; // How many frames to save. != 0
    const float gravity = -9.81;
    const float elasticity = 0.95;
    const int numParticles = 10;
    const float radius = 1; // For smoothing function
    const double timeStep = 1.0/60.0;
    const int mass = 1; // Particle mass


    const box bounds = {100, 0, 100, 0, 0, 0}; //xmax, xmin, ymax, ymin, zmax, zmin

    // ===================== Vector arrays ================================
    vector<double> position(numParticles * numDimensions); // x,y,z interleaved * numParticles
    vector<double> velocity(numParticles * numDimensions); // vx,vy,vz interleaved * numParticles

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
    for (int step = 0; step < maxNumSteps; step++) {
        for (int i = 0; i <numParticles; i++) {
            apply_gravity(i, velocity, timeStep, gravity);
            update_particle_position(i, position, velocity, timeStep);
            bounce_particle(i, position, velocity, bounds, elasticity, radius);
        }

        // Save every saveFrames
        if (step % saveFrames == 0) {
            write_particle_positions_to_csv(csvFile, position, step);
            write_frame_to_gsd(handle, step, numParticles, bounds, position, velocity);
        }
    }
    
    csvFile.close();
    gsd_close(&handle);
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

// Smoothing kernel for calculating density
// Uses a simple (r-x)^3 kernel
float smoothing_kernel(float dist, float radius)
{
    if (dist >= radius) return 0;

    float volume = (M_PI * pow(radius, 4) / 6);
    float influence = radius - dist;
    return influence / volume;
}

float smooth_kernal_derivative(float dist, float radius)
{
    if (dist >= radius) return 0;

    float scale = 12 / (pow(radius, 4) * M_PI);
    float influence = dist - radius;
    return influence * scale;
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

// For vectors
int get_z_index(int i)
{
    return 0;
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
