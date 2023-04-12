//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

void Fluid::transferVelocitiesFromGrid() {
    // Transfer velocities from the grid faces to neighboring particles
    // First, transfer X velocities from the grid to particles
    float yShift = 0.5 * spacing;
    float zShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor(x / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor((y - yShift) / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor((z - zShift) / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4, w5, w6, w7, w8)
        float deltaX = x - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        float deltaZ = (z - zShift) - cellZCoordinate * spacing;
        // Now construct the corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
        float w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
        // From here, compute the cell corners responsible for each weight
        int corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
        int corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
        int corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        // Proceed with transfer from the grid to particle
        int offset = gridHeight * gridWidth;
        // Determine whether contributing corners on grid are valid or not
        float valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
        float valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
        float valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
        float valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
        float valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
        float valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
        float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
        float velocity = particles[i].vx;
        if (d > 0.0)    {
            float picVelocity = (valid1 * w1 * xVelocities[corner1] + valid2 * w2 * xVelocities[corner2] + valid3 * w3 * xVelocities[corner3] + valid4 * w4 * xVelocities[corner4] + valid5 * w5 * xVelocities[corner5] + valid6 * w6 * xVelocities[corner6] + valid7 * w7 * xVelocities[corner7] + valid8 * w8 * xVelocities[corner8]) / d;
            float weightedVelocityChanges = (valid1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + valid2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + valid3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + valid4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4]) + valid5 * w5 * (xVelocities[corner5] - prevXVelocities[corner5]) + valid6 * w6 * (xVelocities[corner6] - prevXVelocities[corner6]) + valid7 * w7 * (xVelocities[corner7] - prevXVelocities[corner7]) + valid8 * w8 * (xVelocities[corner8] - prevXVelocities[corner8])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particles[i].vx = (PIC * picVelocity + (1 - PIC) * flipVelocity);
        }
    }

    // Next, transfer Y velocities from the grid to particles
    float xShift = 0.5 * spacing;
    zShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor((x - xShift) / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor(y / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor((z - zShift) / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4, w5, w6, w7, w8)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = y - cellYCoordinate * spacing;
        float deltaZ = (z - zShift) - cellZCoordinate * spacing;
        // Now construct the corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
        float w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
        // From here, compute the cell corners responsible for each weight
        int corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
        int corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
        int corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 1));
        int corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        // Proceed with transfer from the grid to particle
        int offset = gridWidth;
        // Determine whether contributing corners on grid are valid or not
        float valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
        float valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
        float valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
        float valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
        float valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
        float valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
        float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
        float velocity = particles[i].vy;
        if (d > 0.0)    {
            float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4] + valid5 * w5 * yVelocities[corner5] + valid6 * w6 * yVelocities[corner6] + valid7 * w7 * yVelocities[corner7] + valid8 * w8 * yVelocities[corner8]) / d;
            float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4]) + valid5 * w5 * (yVelocities[corner5] - prevYVelocities[corner5]) + valid6 * w6 * (yVelocities[corner6] - prevYVelocities[corner6]) + valid7 * w7 * (yVelocities[corner7] - prevYVelocities[corner7]) + valid8 * w8 * (yVelocities[corner8] - prevYVelocities[corner8])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particles[i].vy = (PIC * picVelocity + (1 - PIC) * flipVelocity);
        }
    }

    // Next, transfer Z velocities from the grid to particles
    xShift = 0.5 * spacing;
    yShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor((x - xShift) / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor((y - yShift) / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor(z / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4, w5, w6, w7, w8)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        float deltaZ = (z) - cellZCoordinate * spacing;
        // Now construct the corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1 - (deltaZ / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1 - (deltaZ / spacing));
        float w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
        float w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
        float w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
        // From here, compute the cell corners responsible for each weight
        int corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
        int corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
        int corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
        int corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        int corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
        // Proceed with transfer from the grid to particle
        int offset = 1;
        // Determine whether contributing corners on grid are valid or not
        float valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
        float valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
        float valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
        float valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
        float valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
        float valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
        float valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
        float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
        float velocity = particles[i].vz;
        if (d > 0.0)    {
            float picVelocity = (valid1 * w1 * zVelocities[corner1] + valid2 * w2 * zVelocities[corner2] + valid3 * w3 * zVelocities[corner3] + valid4 * w4 * zVelocities[corner4] + valid5 * w5 * zVelocities[corner5] + valid6 * w6 * zVelocities[corner6] + valid7 * w7 * zVelocities[corner7] + valid8 * w8 * zVelocities[corner8]) / d;
            float weightedVelocityChanges = (valid1 * w1 * (zVelocities[corner1] - prevZVelocities[corner1]) + valid2 * w2 * (zVelocities[corner2] - prevZVelocities[corner2]) + valid3 * w3 * (zVelocities[corner3] - prevZVelocities[corner3]) + valid4 * w4 * (zVelocities[corner4] - prevZVelocities[corner4]) + valid5 * w5 * (zVelocities[corner5] - prevZVelocities[corner5]) + valid6 * w6 * (zVelocities[corner6] - prevZVelocities[corner6]) + valid7 * w7 * (zVelocities[corner7] - prevZVelocities[corner7]) + valid8 * w8 * (zVelocities[corner8] - prevZVelocities[corner8])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particles[i].vz = (PIC * picVelocity + (1 - PIC) * flipVelocity);
        }
    }
}