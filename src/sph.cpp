// Created by Markus Zuegg
// Main cpp script to run sph simulation
// Main simulation loop will be housed here
#include <bits/stdc++.h>

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

// ========================== Functions ========================================

int main(int argc, char** argv)
{
    const float gravity = 9.81;
    const float elasticity = 0.95;
    const int numParticles = 5;
    const float radius = 1;

    vector<double> position(numParticles * 2); // Consists of x,y * numParticles
    vector<double> velocity(numParticles * 2); // Consists of x,y * numParticles
    
    
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
        velocity[x] = -velocity[x] * elasticity;
    }
    if (position[x] + radius > bounds.xmax) {
        position[x] = bounds.xmax - radius;
        velocity[x] = -velocity[x] * elasticity;
    }
    // Check walls in y
    if (position[y] - radius < bounds.ymin) {
        position[y] = bounds.ymin + radius;
        velocity[y] = -velocity[y] * elasticity;
    }
    if (position[y] + radius > bounds.ymax) {
        position[y] = bounds.ymax - radius;
        velocity[y] = -velocity[y] * elasticity;
    }
}
