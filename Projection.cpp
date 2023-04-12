//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

void Fluid::projectGS() {
    // Gauss-Seidel solver (without vectorization)
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    prevZVelocities = zVelocities;
    int numIterations = 200;
    float overRelaxation = 1.9;
    for (int i = 0; i < numIterations; i++) {
        // Note that we go from 1 to gridDimension - 2 to make sure there are always valid velocity components to pull from
        for (int cellXCoordinate = 1; cellXCoordinate < gridLength - 1; cellXCoordinate++)    {
            for (int cellYCoordinate = 1; cellYCoordinate < gridHeight - 1; cellYCoordinate++)    {
                for (int cellZCoordinate = 1; cellZCoordinate < gridWidth - 1; cellZCoordinate++)   {
                    if (cellType[IX(cellXCoordinate, cellYCoordinate, cellZCoordinate)] == FLUID) {
                        // Define top, bottom, left and right cell coordinates
                        int cellCoordinate = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
                        int leftCellCoordinate = IX(cellXCoordinate - 1, cellYCoordinate, cellZCoordinate);
                        int rightCellCoordinate = IX(cellXCoordinate + 1, cellYCoordinate, cellZCoordinate);
                        int bottomCellCoordinate = IX(cellXCoordinate, cellYCoordinate - 1, cellZCoordinate);
                        int topCellCoordinate = IX(cellXCoordinate, cellYCoordinate + 1, cellZCoordinate);
                        int backCellCoordinate = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate - 1);
                        int frontCellCoordinate = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate + 1);
                        // Extract velocity components
                        float leftXVelocity = xVelocities[cellCoordinate];
                        float rightXVelocity = xVelocities[rightCellCoordinate];
                        float bottomYVelocity = yVelocities[cellCoordinate];
                        float topYVelocity = yVelocities[topCellCoordinate];
                        float backZVelocity = zVelocities[cellCoordinate];
                        float frontZVelocity = zVelocities[frontCellCoordinate];
                        // Initialize divergence
                        float divergence = 0.0;
                        // Compute multiplication factor (inspired by Ten Minute Physics' explanation to account for solid cells)
                        float leftMultiplicationFactor =
                                cellType[leftCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float rightMultiplicationFactor =
                                cellType[rightCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float topMultiplicationFactor =
                                cellType[topCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float bottomMultiplicationFactor =
                                cellType[bottomCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float backMultiplicationFactor =
                                cellType[backCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float frontMultiplicationFactor =
                                cellType[frontCellCoordinate] != SOLID ? 1.0 : 0.0;
                        float normalizationFactor =
                                leftMultiplicationFactor + rightMultiplicationFactor + topMultiplicationFactor +
                                bottomMultiplicationFactor + backMultiplicationFactor + frontMultiplicationFactor;
                        if (normalizationFactor == 0) continue;
                        // Finally, compute divergence and modify velocities considering all multiplication factors
                        divergence = (rightXVelocity - leftXVelocity + topYVelocity - bottomYVelocity + frontZVelocity - backZVelocity);
                        // Account for the particle density by subtracting the compression factor
                        float stiffnessCoefficient = VOLUME_PRESERVATION_TERM;
                        if (initialParticleDensity > 0) {
                            float compression = cellParticleDensity[cellCoordinate] - initialParticleDensity;
                            //if (compression > 0) divergence -= (stiffnessCoefficient * compression) / spacing;
                        }
                        divergence *= -1 / normalizationFactor;
                        xVelocities[cellCoordinate] -= divergence * leftMultiplicationFactor * overRelaxation;
                        xVelocities[rightCellCoordinate] += divergence * rightMultiplicationFactor * overRelaxation;
                        yVelocities[cellCoordinate] -= divergence * bottomMultiplicationFactor * overRelaxation;
                        yVelocities[topCellCoordinate] += divergence * topMultiplicationFactor * overRelaxation;
                        zVelocities[cellCoordinate] -= divergence * backMultiplicationFactor * overRelaxation;
                        zVelocities[frontCellCoordinate] += divergence * frontMultiplicationFactor * overRelaxation;
                    }
                }
            }
        }
    }
}

void Fluid::computeCellDensities() {
    // First, reset the current cellParticleDensity grid
    std::fill(cellParticleDensity.begin(), cellParticleDensity.end(), 0.0);
    // Since density is stored at cell centers, we need to shift both coordinates by half the grid spacing
    float xShift = 0.5f * spacing;
    float yShift = 0.5f * spacing;
    float zShift = 0.5f * spacing;
    // Now, iterate over all particles
    for (int i = 0; i < numParticles; i++) {
        // for each particle, compute the cell coordinate and contribution of weight through interpolation
        int cellXCoordinate = Utils::clamp(floor((particles[i].x - xShift) / spacing), 0, gridLength - 1);
        int cellYCoordinate = Utils::clamp(floor((particles[i].y - yShift) / spacing), 0, gridHeight - 1);
        int cellZCoordinate = Utils::clamp(floor((particles[i].z - zShift) / spacing), 0, gridWidth - 1);
        float deltaX = (particles[i].x - xShift) - spacing * cellXCoordinate;
        float deltaY = (particles[i].y - yShift) - spacing * cellYCoordinate;
        float deltaZ = (particles[i].z - zShift) - spacing * cellZCoordinate;
        float w1 = (1 - (deltaX / spacing)) * (1 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w2 = (deltaX / spacing) * (1 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w4 = (1 - (deltaX / spacing)) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w5 = (1 - (deltaX / spacing)) * (1 - (deltaY / spacing)) * (deltaZ / spacing);
        float w6 = (deltaX / spacing) * (1 - (deltaY / spacing)) * (deltaZ / spacing);
        float w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
        float w8 = (1 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
        int corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
        int corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
        int corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        cellParticleDensity[corner1] += w1;
        cellParticleDensity[corner2] += w2;
        cellParticleDensity[corner3] += w3;
        cellParticleDensity[corner4] += w4;
        cellParticleDensity[corner5] += w5;
        cellParticleDensity[corner6] += w6;
        cellParticleDensity[corner7] += w7;
        cellParticleDensity[corner8] += w8;
    }

    if (initialParticleDensity == 0)    {
        // compute when the simulation starts
        float accumulator = 0;
        int fluidCellCount = 0;
        for (int i = 0; i < numCells; i++)  {
            if (cellType[i] == FLUID)   {
                fluidCellCount++;
                accumulator += cellParticleDensity[i];
            }
        }
        if (fluidCellCount) initialParticleDensity = accumulator / fluidCellCount;
    }
}