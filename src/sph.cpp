// Created by Markus Zuegg
// Main cpp script to run sph simulation
// Main simulation loop will be housed here
#include <bits/stdc++.h>
#include <string>
#include <random>

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
H5::H5File initialize_hdf5_file(const std::string& filename, const box& bounds);
void save_frame_to_hdf5(H5::H5File& file, const vector<double>& position, int frame_num, const std::string& prefix);

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    argc--;
    argv++;
    if (argc < 0) return -1; // Check for filename arg

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
    char* fileSuffix = "_particle_positions.h5";
    const std::string FILE_NAME(strcat(argv[0], fileSuffix));
    const std::string DATASET_PREFIX("frame_");
    H5::H5File file = initialize_hdf5_file(FILE_NAME, bounds);

    // Simulation loop
    for (int step = 0; step < maxNumSteps; step++) {
        for (int i = 0; i <numParticles; i++) {
            velocity[i*2 + 1] += gravity * timeStep; // vy += g * dt
            position[i*2] += velocity[i*2] * timeStep; // x += vx * dt
            position[i*2 + 1] += velocity[2*i + 1] * timeStep; // y += vy * dt
            bounce_particle(i, position, velocity, bounds, elasticity, radius);
        }
        save_frame_to_hdf5(file, position, step, DATASET_PREFIX);
    }
    
    return 0;
}

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
