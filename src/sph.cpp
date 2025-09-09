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
void apply_gravity(int i, vector<double> &velocity, float timeStep, float gravity);
void calculate_denisty(int i, vector<double> &densities, vector<double> &position,
        float mass, int numParticles, float radius);
float density_to_pressure(double density, float targetDensity, float pressureMulti);
void calculate_pressure(int i, double &force_x, double &force_y, vector<double> &postion,
                        vector<double> &densities, int numParticles, int mass, float radius,
                        float targetDensity, float pressureMulti);

int get_x_index(int i);
int get_y_index(int i);
int get_z_index(int i);
float distance_from_two_points(double x1, double y1, double x2, double y2);

float ploy6_kernal(float dist, float radius);
float spikey_kernal(float dist, float radius);
float spikey_kernal_derivative(float dist, float radius);

void init_random_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float radius);
void init_square_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float spacing);

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

    // ============= Simulation params ==================
    const int maxNumSteps = 1500;
    const int saveFrames = 1; // How many frames to save. != 0
    const int numParticles = 400;
    const double timeStep = 1.0/60.0;
    
    const float gravity = 0;
    const float elasticity = 0.85;
    
    const int mass = 1; // Particle mass
    const float radius = 0.1; // for particles

    const float smoothingRadius = 1.2;
    const float targetDensity = 1.5;
    const float pressureMulti = 1;


    const box bounds = {10, 0, 10, 0, 0, 0}; //xmax, xmin, ymax, ymin, zmax, zmin

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
    float spacing = 2.0;
    init_random_particle_pos(position, velocity, numParticles, bounds, radius); // random starting pos
    for (int step = 0; step < maxNumSteps; step++) {
        for (int i = 0; i < numParticles; i++) {
            // Apply gravity and calc denisties for all particles
            // TODO make it use predictive positions
            apply_gravity(i, velocity, timeStep, gravity);
            calculate_denisty(i, densities, position, mass, numParticles, smoothingRadius);
            
        }

        for (int i = 0; i < numParticles; i++) {
            // Calc pressure
            double forceX = 0.0f;
            double forceY = 0.0f;
            calculate_pressure(i, forceX, forceY, position, densities,
                    numParticles, mass, smoothingRadius, targetDensity, pressureMulti);
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
    force_x = 0.0f;
    force_y = 0.0f;

    // Get current particle position
    double x1 = postion[get_x_index(i)];
    double y1 = postion[get_y_index(i)];

    // Get density and pressure for 1st particle
    double density1 = densities[i];
    float pressure1 = density_to_pressure(density1, targetDensity, pressureMulti);
    if (density1 <= 0.0f) return;

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

void init_square_particle_pos(vector<double> &position, vector<double> &velocity,
        int numParticles, const box &bounds, float spacing)
{
    float centerX = 0.5 * (bounds.xmax + bounds.xmin);
    float centerY = 0.5 * (bounds.ymax + bounds.ymin);

    int gridSize = static_cast<int>(ceil(sqrt(numParticles)));
    float halfSpan = 0.5 * (gridSize - 1) * spacing;

    int count = 0;
    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            if (count >= numParticles) break;

            int x = get_x_index(count);
            int y = get_y_index(count);

            position[x] = centerX - halfSpan + col * spacing;
            position[y] = centerY - halfSpan + row * spacing;

            velocity[x] = 0.0;
            velocity[y] = 0.0;

            count++;
        }
    }
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
