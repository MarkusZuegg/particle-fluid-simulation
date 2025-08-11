// Created by Markus Zuegg
// Main cpp script to run sph simulation
// Main simulation loop will be housed here
#include <bits/stdc++.h>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <H5Cpp.h>

using std::vector;

// =========================== Const defs ======================================

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
        float timeStep, float gravity);

std::ofstream open_csv_file(const std::string& filename, int numParticles, const box& bounds);
void write_particle_positions_to_csv(std::ofstream &file, const vector<double>& positions);

H5::H5File initialize_hdf5_file(const std::string& filename, const box& bounds);
void save_frame_to_hdf5(H5::H5File& file, const vector<double>& position, int frame_num,
        const std::string& prefix);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 1) return -1; // Check for filename arg

    // ============= simulation params ==================
    const int maxNumSteps = 500;
    const float gravity = -9.81;
    const float elasticity = 0.95;
    const int numParticles = 5;
    const float radius = 1;
    const double timeStep = 1/60;

    vector<double> position(numParticles * 2); // x,y interleaved * numParticles
    vector<double> velocity(numParticles * 2); // vx,vy interleaved * numParticles
    
    const box bounds = {100, 0, 100, 0}; //xmax, xmin, ymax, ymin

    // ============== HDF5 setup ========================
    // char fileSuffix[] = "_particle_positions.h5";
    // char* jobID = strdup(argv[0]);
    
    std::string fileName = std::string(argv[0]) + ".csv";
    std::ofstream csvFile = open_csv_file(fileName, numParticles, bounds);

    // Simulation loop
    for (int step = 0; step < maxNumSteps; step++) {
        for (int i = 0; i <numParticles; i++) {
            update_particle_position(i, position, velocity, timeStep, gravity);
            bounce_particle(i, position, velocity, bounds, elasticity, radius);
        }

        write_particle_positions_to_csv(csvFile, position);
    }
    
    csvFile.close();
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

void update_particle_position(int i, vector<double> &position, vector<double> &velocity,
        float timeStep, float gravity)
{
    // Get particle index for x and y
    int x = 2 * i;
    int y = 2 * i + 1;

    velocity[y] += gravity * timeStep; // vy += g * dt
    position[x] += velocity[x] * timeStep; // x += vx * dt
    position[y] += velocity[y] * timeStep; // y += vy * dt
}

//=========================== csv saving =======================================

// 
std::ofstream open_csv_file(const std::string& filename, int numParticles, const box& bounds)
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << numParticles << "\n";
        file << bounds.xmin << "\n";
        file << bounds.xmax << "\n";
        file << bounds.ymin << "\n";
        file << bounds.ymax << "\n";
        return file;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void write_particle_positions_to_csv(std::ofstream &file, const vector<double>& positions)
{
    if (file.is_open()) {
        for (size_t i = 0; i < positions.size(); ++i) {
            file << positions[i] << (i == positions.size() - 1 ? "" : ",");
        }
        file << "\n";
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

//========================== hdf5 stuff ========================================

H5::H5File initialize_hdf5_file(const std::string& filename, const box& bounds)
{
    H5::H5File file(filename, H5F_ACC_TRUNC);
    double bounds_data[4] = {bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax};

    hsize_t dims[1] = {4};
    H5::DataSpace attr_dataspace(1, dims);

    H5::Attribute attribute = file.createAttribute("bounds", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    attribute.write(H5::PredType::NATIVE_DOUBLE, bounds_data);

    std::cout << "HDF5 file created and bounds attribute written." << std::endl;
    return file;
}

void save_frame_to_hdf5(H5::H5File& file, const vector<double>& position, int frame_num, const std::string& prefix)
{
    hsize_t dims[1] = {position.size()};
    H5::DataSpace dataspace(1, dims);
    std::string dataset_name = prefix + std::to_string(frame_num);
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(position.data(), H5::PredType::NATIVE_DOUBLE);
}
