//
// Created by Anchit Mishra on 2023-04-08.
//

#include "Fluid.h"

using std::tuple;

Fluid::Fluid() {
    // Part 1: Initialize grid stuff
    spacing = fmax(fmax(WIDTH / gridWidth, LENGTH / gridLength), HEIGHT / gridHeight);
    // Initialize the grid velocities
    xVelocities.resize(numCells, 0.0);
    yVelocities.resize(numCells, 0.0);
    zVelocities.resize(numCells, 0.0);
    // Initialize the marker array
    xMarker.resize(numCells, INT_MAX);
    yMarker.resize(numCells, INT_MAX);
    zMarker.resize(numCells, INT_MAX);
    // Initialize the previous grid velocities
    prevXVelocities.resize(numCells, 0.0);
    prevYVelocities.resize(numCells, 0.0);
    prevZVelocities.resize(numCells, 0.0);
    // Initialize the delta grid velocities
    xR.resize(numCells, 0.0);
    yR.resize(numCells, 0.0);
    zR.resize(numCells, 0.0);
    // Initialize the cell type matrix
    cellType.resize(numCells, 0);
    // Initialize the cell particle density matrix
    cellParticleDensity.resize(numCells, 0.0);
    initialParticleDensity = 0;
    // Initialize the unbounded spatial hash table
    spatialHashTableSize = spatialHashGridX * spatialHashGridY * spatialHashGridZ;
    spatialHashTable.resize(spatialHashTableSize + 1, 0);
    spatialHashParticles.resize(numParticles);

    // Part 2: Initialize particle positions and velocities
    particles.resize(numParticles);
    for (int i = 0; i < numParticles; i++)  {
        particles[i].x = 0;
        particles[i].y = 0;
        particles[i].z = 0;
        particles[i].vx = 0;
        particles[i].vy = 0;
        particles[i].vz = 0;
    }



    // Now that basic initialization is done, position the particles in a starting configuration for dam-break
    int index = 0;
    for (int i = 0; i < particleMassLength; i++)    {
        for (int j = 0; j < particleMassHeight; j++)    {
            for (int k = 0; k < particleMassWidth; k++) {
                particles[index].x = 0.2 * LENGTH + SPACING + particleRadius + 2 * particleRadius * i + (j % 2 == 0 ? 0 : particleRadius);
                particles[index].y = 0.05 * HEIGHT + SPACING + particleRadius + 1.73 * particleRadius * j;
                particles[index++].z = 0.05 * WIDTH + SPACING + particleRadius + 2 * particleRadius * k;
            }
        }
    }

    containerWallXVelocity = 0;
    containerWallYVelocity = 0;
    containerWallZVelocity = 0;

}

int Fluid::IX(int x, int y, int z) {
    return x * gridHeight * gridWidth + y * gridWidth + z;
}

void Fluid::simulateFluid() {
    // Implement simulation cycle here - transfer to grid, solve incompressibility, transfer from grid.
    //std::cout << "Simulating fluid...\n";
    //std::cout << "Advecting...\n";
    advect();
    //std::cout << "Penalties...\n";
    if (PENALTY) {
        spatialHashing();
        detectParticleCollisions();
    }
    //std::cout << "Detecting boundaries...\n";
    detectBoundaryCollisions();
    //std::cout << "TTG...\n";
    markFluidCells();
    transferVelocitiesToGrid();
    //std::cout << "Extrapolating...\n";
    extrapolateVelocities();
    //std::cout << "Projecting...\n";
    if (BETTER_PROJECTION)  {
        projectPCG();
    } else {
        computeCellDensities();
        projectGS();
    }
    //std::cout << "Extrapolating...\n";
    extrapolateVelocities();
    //std::cout << "TTP...\n";
    transferVelocitiesFromGrid();
    //std::cout << "Simulation step complete.\n";
}
