//
// Created by Anchit Mishra on 2023-04-08.
//

#ifndef FLUIDSIM_3D_FLUID_H
#define FLUIDSIM_3D_FLUID_H

#include <iostream>
#include <cmath>
#include <vector>
#include "Utils.h"
#include "Particle.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

#define LENGTH 1280
#define HEIGHT 640
#define WIDTH 640
#define SPACING 10.0
#define PARTICLE_SIZE 0.25
#define DENSITY 1000.0
#define SOLID 2
#define FLUID 0
#define EMPTY 1
#define PUSH_PENALTY 0.1
#define VOLUME_PRESERVATION_TERM 1
#define BETTER_PROJECTION false
#define BETTER_INTEGRATION true
#define TIMESTEP 1.0/60.0
#define SUBSTEPS 2
#define PENALTY false
#define PIC 0.1

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;

class Fluid {
public:
    // Fluid should contain a grid of specified dimensions
    float spacing; // This is the length of one cell in the grid
    int gridLength = floor(LENGTH / SPACING) + 1;
    int gridHeight = floor(HEIGHT / SPACING) + 1;
    int gridWidth = floor(WIDTH / SPACING) + 1;
    int numCells = gridLength * gridHeight * gridWidth;

    // The grid is MAC - it stores velocities on the cube faces, pressure in the center
    // Current velocities
    std::vector<float> xVelocities;
    std::vector<float> yVelocities;
    std::vector<float> zVelocities;
    // Marker array for extrapolation
    std::vector<int> xMarker;
    std::vector<int> yMarker;
    std::vector<int> zMarker;
    // Previous velocities
    std::vector<float> prevXVelocities;
    std::vector<float> prevYVelocities;
    std::vector<float> prevZVelocities;
    // Normalization coefficients (for interpolation)
    std::vector<float> xR;
    std::vector<float> yR;
    std::vector<float> zR;
    // Each grid cell should also be marked as being solid, fluid or empty
    std::vector<int> cellType;
    std::vector<float> cellParticleDensity;
    // We also need to track the initial particle density;
    float initialParticleDensity = 0.0; // We want roughly 4 particles per cell (2 along each dimension)
    // Spatial hashing table (unbounded grid)
    std::vector<int> spatialHashTable;
    std::vector<int> spatialHashParticles;
    int spatialHashTableSize;

    // Next up, we store penalty grid-related information for the fluid
    float particleRadius = SPACING * PARTICLE_SIZE;
    float spatialHashGridSpacing = 2.2 * particleRadius;
    int spatialHashGridX = floor(LENGTH / spatialHashGridSpacing) + 1;
    int spatialHashGridY = floor(HEIGHT / spatialHashGridSpacing) + 1;
    int spatialHashGridZ = floor(WIDTH / spatialHashGridSpacing) + 1;
    int particleMassLength = floor((0.3 * (LENGTH) - (2.0 * particleRadius) - (2.0 * SPACING)) / (2.0 * particleRadius));
    int particleMassHeight = floor((0.9 * (HEIGHT) - (2.0 * particleRadius) - (2.0 * SPACING)) / (2.0 * particleRadius));
    int particleMassWidth = floor((0.9 * (WIDTH) - (2.0 * particleRadius) - (2.0 * SPACING)) / (2.0 * particleRadius));
    int numParticles = particleMassHeight * particleMassLength * particleMassWidth;
    // store particles
    std::vector<Particle> particles;
    // Faking container wall velocities for now
    float containerWallXVelocity;
    float containerWallYVelocity;
    float containerWallZVelocity;

    // Finally, we come to functions
    Fluid();
    void simulateFluid();
    int IX(int x, int y, int z);
    // Advection
    void advect();
    void detectBoundaryCollisions();
    void extrapolateVelocities();
    // Transfer
    void markFluidCells();
    void transferVelocitiesToGrid();
    void transferVelocitiesFromGrid();
    // Projection
    void projectGS();
    void projectPCG();
    void computeCellDensities();
    // Spatial Hashing (for collision detection)
    int spatialHashFunction(int x, int y, int z);
    void spatialHashing();
    void detectParticleCollisions();
    // Surface reconstruction
    float distance(std::tuple<int, int, int> x, std::tuple<float, float, float> p);
};


#endif //FLUIDSIM_3D_FLUID_H
