//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

// Hash function used as per description in 10 minute physics tutorial on spatial hashing
int Fluid::spatialHashFunction(int x, int y, int z) {
    return abs((x * 92837111) ^ (y * 689287499) ^ (z * 283923481));
}

void Fluid::spatialHashing() {
    std::fill(spatialHashTable.begin(), spatialHashTable.end(), 0);
    std::fill(spatialHashParticles.begin(), spatialHashParticles.end(), -1);
    // Compute which cell each particle falls into
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = Utils::clamp(floor(particles[i].x / spatialHashGridSpacing), 0, spatialHashGridX - 1);
        int cellYCoordinate = Utils::clamp(floor(particles[i].y / spatialHashGridSpacing), 0, spatialHashGridY - 1);
        int cellZCoordinate = Utils::clamp(floor(particles[i].z / spatialHashGridSpacing), 0, spatialHashGridZ - 1);
        int cellIndex = cellXCoordinate * spatialHashGridY * spatialHashGridZ + cellYCoordinate * spatialHashGridZ + cellZCoordinate;
        // increase particle count at the specified index
        spatialHashTable[cellIndex]++;
    }

    // Iterate over the table and compute partial sums
    int accumulator = 0;
    for (int i = 0; i < spatialHashTableSize; i++)  {
        accumulator += spatialHashTable[i];
        spatialHashTable[i] = accumulator;
    }
    // final entry stores the last accumulator value
    spatialHashTable[spatialHashTableSize] = accumulator;

    // Now, populate the particle vector corresponding to the spatial hash
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = Utils::clamp(floor(particles[i].x / spatialHashGridSpacing), 0, spatialHashGridX - 1);
        int cellYCoordinate = Utils::clamp(floor(particles[i].y / spatialHashGridSpacing), 0, spatialHashGridY - 1);
        int cellZCoordinate = Utils::clamp(floor(particles[i].z / spatialHashGridSpacing), 0, spatialHashGridZ - 1);
        int cellIndex = cellXCoordinate * spatialHashGridY * spatialHashGridZ + cellYCoordinate * spatialHashGridZ + cellZCoordinate;
        spatialHashTable[cellIndex]--;
        // store the index of the particle
        spatialHashParticles[spatialHashTable[cellIndex]] = i;
    }

}

void Fluid::detectParticleCollisions() {
    // of course, iterate over all particles
    for (int iter = 0; iter < 2; iter++) {
        for (int i = 0; i < numParticles; i++) {
            // we want to detect particles within a distance of 2 * particleRadius
            int cellXCoordinate = Utils::clamp(floor(particles[i].x / spatialHashGridSpacing), 0, spatialHashGridX - 1);
            int cellYCoordinate = Utils::clamp(floor(particles[i].y / spatialHashGridSpacing), 0, spatialHashGridY - 1);
            int cellZCoordinate = Utils::clamp(floor(particles[i].z / spatialHashGridSpacing), 0, spatialHashGridZ - 1);
            int lowX = fmax(cellXCoordinate - 1, 0);
            int lowY = fmax(cellYCoordinate - 1, 0);
            int lowZ = fmax(cellZCoordinate - 1, 0);
            int highX = fmin(cellXCoordinate + 1, spatialHashGridX - 1);
            int highY = fmin(cellYCoordinate + 1, spatialHashGridY - 1);
            int highZ = fmin(cellZCoordinate + 1, spatialHashGridZ - 1);
            for (int j = lowX; j < highX + 1; j++) {
                for (int k = lowY; k < highY + 1; k++) {
                    for (int l = lowZ; l < highZ + 1; l++) {
                        int tableIndex = j * spatialHashGridY * spatialHashGridZ + k * spatialHashGridZ + l;
                        // check if particles are present in this index of the table
                        int numParticlesInCell = spatialHashTable[tableIndex + 1] - spatialHashTable[tableIndex];
                        if (numParticlesInCell > 0) {
                            // check particles
                            for (int a = spatialHashTable[tableIndex];
                                 a < spatialHashTable[tableIndex] + numParticlesInCell; a++) {
                                int particle = spatialHashParticles[a];
                                if (particle != i) {
                                    // compute distance
                                    float x = (particles[i].x - particles[particle].x);
                                    float y = (particles[i].y - particles[particle].y);
                                    float z = (particles[i].z - particles[particle].z);
                                    float distance2 = (pow(x, 2) + pow(y, 2) + pow(z, 2));
                                    float distance = sqrt(distance2);
                                    if (distance2 < 4.0 * particleRadius * particleRadius && distance2 > 0) {
                                        // apply penalties, push particles (spring-like penalties)
                                        float penaltyFactor = 0.5 * (2.0 * particleRadius - distance) / distance;
                                        particles[i].x += penaltyFactor * x;
                                        particles[i].y += penaltyFactor * y;
                                        particles[i].z += penaltyFactor * z;
                                        particles[particle].x -= penaltyFactor * x;
                                        particles[particle].y -= penaltyFactor * y;
                                        particles[particle].z -= penaltyFactor * z;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}