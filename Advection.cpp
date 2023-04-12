//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

void Fluid::detectBoundaryCollisions() {
    // First, set a minimum and maximum limit for particles to be located at
    float minX = spacing + particleRadius;
    float minY = spacing + particleRadius;
    float minZ = spacing + particleRadius;
    float maxX = spacing * (gridLength - 1) - particleRadius;
    float maxY = spacing * (gridHeight - 1) - particleRadius;
    float maxZ = spacing * (gridWidth - 1) - particleRadius;
    // Iterate over all particles, see whether they are colliding with any boundaries
    for (int i = 0; i < numParticles; i++)  {
        // Clamp the X and Y coordinates of the particle between boundary coordinates
        if (particles[i].x < minX || particles[i].x > maxX) particles[i].vx = (containerWallXVelocity);
        if (particles[i].y < minY || particles[i].y > maxY) particles[i].vy = (containerWallYVelocity);
        if (particles[i].z < minZ || particles[i].z > maxZ) particles[i].vz = (containerWallZVelocity);
        particles[i].x = (Utils::clamp(particles[i].x, minX, maxX));
        particles[i].y = (Utils::clamp(particles[i].y, minY, maxY));
        particles[i].z = (Utils::clamp(particles[i].z, minZ, maxZ));
    }
}

void Fluid::advect()    {
    float dt = TIMESTEP/SUBSTEPS;
    if (BETTER_INTEGRATION) {
        for (int i = 0; i < numParticles; i++)  {
            // Extract stage 1 values, k1_i
            float k1_x = particles[i].vx;
            float k1_y = particles[i].vy;
            float k1_z = particles[i].vz;
            // Now comes the hard part - derive stage 2 velocities
            float k2_x = 0;
            float k2_y = 0;
            float k2_z = 0;
            float secondXPosition = particles[i].x + 0.5 * dt * k1_x;
            float secondYPosition = particles[i].y + 0.5 * dt * k1_y;
            float secondZPosition = particles[i].z + 0.5 * dt * k1_z;
            float yShift = 0.5 * spacing;
            float zShift = 0.5 * spacing;
            // Interpolate the velocity field to get velocities at this position
            int cellXCoordinate = floor(secondXPosition / spacing);
            int cellYCoordinate = floor((secondYPosition - yShift) / spacing);
            int cellZCoordinate = floor((secondZPosition - zShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            float deltaX = secondXPosition - cellXCoordinate * spacing;
            float deltaY = (secondYPosition - yShift) - cellYCoordinate * spacing;
            float deltaZ = (secondZPosition - zShift) - cellZCoordinate * spacing;
            // Now construct the corner weights
            float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            float w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
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
                k2_x = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            float xShift = 0.5 * spacing;
            zShift = 0.5 * spacing;
            // Do the same for stage 2 y velocities
            cellXCoordinate = floor((secondXPosition - xShift) / spacing);
            cellYCoordinate = floor(secondYPosition / spacing);
            cellZCoordinate = floor((secondZPosition - zShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (secondXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = secondYPosition - cellYCoordinate * spacing;
            deltaZ = (secondZPosition - zShift) - cellZCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
            w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
            corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
            corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            // Proceed with transfer from the grid to particle
            offset = gridWidth;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
            valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
            valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
            valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
            valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
            valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
            velocity = particles[i].vy;
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4] + valid5 * w5 * yVelocities[corner5] + valid6 * w6 * yVelocities[corner6] + valid7 * w7 * yVelocities[corner7] + valid8 * w8 * yVelocities[corner8]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4]) + valid5 * w5 * (yVelocities[corner5] - prevYVelocities[corner5]) + valid6 * w6 * (yVelocities[corner6] - prevYVelocities[corner6]) + valid7 * w7 * (yVelocities[corner7] - prevYVelocities[corner7]) + valid8 * w8 * (yVelocities[corner8] - prevYVelocities[corner8])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k2_y = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            xShift = 0.5 * spacing;
            yShift = 0.5 * spacing;
            // Do the same for stage 2 z velocities
            cellXCoordinate = floor((secondXPosition - xShift) / spacing);
            cellYCoordinate = floor((secondYPosition - yShift) / spacing);
            cellZCoordinate = floor(secondZPosition / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (secondXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = (secondYPosition - yShift) - cellYCoordinate * spacing;
            deltaZ = secondZPosition - cellZCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
            w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
            corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
            corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            // Proceed with transfer from the grid to particle
            offset = 1;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
            valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
            valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
            valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
            valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
            valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
            velocity = particles[i].vz;
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * zVelocities[corner1] + valid2 * w2 * zVelocities[corner2] + valid3 * w3 * zVelocities[corner3] + valid4 * w4 * zVelocities[corner4] + valid5 * w5 * zVelocities[corner5] + valid6 * w6 * zVelocities[corner6] + valid7 * w7 * zVelocities[corner7] + valid8 * w8 * zVelocities[corner8]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (zVelocities[corner1] - prevZVelocities[corner1]) + valid2 * w2 * (zVelocities[corner2] - prevZVelocities[corner2]) + valid3 * w3 * (zVelocities[corner3] - prevZVelocities[corner3]) + valid4 * w4 * (zVelocities[corner4] - prevZVelocities[corner4]) + valid5 * w5 * (zVelocities[corner5] - prevZVelocities[corner5]) + valid6 * w6 * (zVelocities[corner6] - prevZVelocities[corner6]) + valid7 * w7 * (zVelocities[corner7] - prevZVelocities[corner7]) + valid8 * w8 * (zVelocities[corner8] - prevZVelocities[corner8])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k2_z = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            // Finally, derive the 3rd stage velocities
            float k3_x = 0;
            float k3_y = 0;
            float k3_z = 0;
            float thirdXPosition = particles[i].x + 0.75 * dt * k2_x;
            float thirdYPosition = particles[i].y + 0.75 * dt * k2_y;
            float thirdZPosition = particles[i].z + 0.75 * dt * k3_z;

            yShift = 0.5 * spacing;
            zShift = 0.5 * spacing;
            // Interpolate the velocity field to get velocities at this position
            cellXCoordinate = floor(thirdXPosition / spacing);
            cellYCoordinate = floor((thirdYPosition - yShift) / spacing);
            cellZCoordinate = floor((thirdZPosition - zShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = thirdXPosition - cellXCoordinate * spacing;
            deltaY = (thirdYPosition - yShift) - cellYCoordinate * spacing;
            deltaZ = (thirdZPosition - zShift) - cellZCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
            w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
            corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
            corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            // Proceed with transfer from the grid to particle
            offset = gridHeight * gridWidth;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
            valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
            valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
            valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
            valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
            valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
            velocity = particles[i].vx;
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * xVelocities[corner1] + valid2 * w2 * xVelocities[corner2] + valid3 * w3 * xVelocities[corner3] + valid4 * w4 * xVelocities[corner4] + valid5 * w5 * xVelocities[corner5] + valid6 * w6 * xVelocities[corner6] + valid7 * w7 * xVelocities[corner7] + valid8 * w8 * xVelocities[corner8]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + valid2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + valid3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + valid4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4]) + valid5 * w5 * (xVelocities[corner5] - prevXVelocities[corner5]) + valid6 * w6 * (xVelocities[corner6] - prevXVelocities[corner6]) + valid7 * w7 * (xVelocities[corner7] - prevXVelocities[corner7]) + valid8 * w8 * (xVelocities[corner8] - prevXVelocities[corner8])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k3_x = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            xShift = 0.5 * spacing;
            zShift = 0.5 * spacing;
            // Do the same for stage 3 y velocities
            cellXCoordinate = floor((thirdXPosition - xShift) / spacing);
            cellYCoordinate = floor(thirdYPosition / spacing);
            cellZCoordinate = floor((thirdZPosition - zShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (thirdXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = thirdYPosition - cellYCoordinate * spacing;
            deltaZ = (thirdZPosition - zShift) - cellZCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
            w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
            corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
            corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            // Proceed with transfer from the grid to particle
            offset = gridWidth;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
            valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
            valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
            valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
            valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
            valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
            velocity = particles[i].vy;
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4] + valid5 * w5 * yVelocities[corner5] + valid6 * w6 * yVelocities[corner6] + valid7 * w7 * yVelocities[corner7] + valid8 * w8 * yVelocities[corner8]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4]) + valid5 * w5 * (yVelocities[corner5] - prevYVelocities[corner5]) + valid6 * w6 * (yVelocities[corner6] - prevYVelocities[corner6]) + valid7 * w7 * (yVelocities[corner7] - prevYVelocities[corner7]) + valid8 * w8 * (yVelocities[corner8] - prevYVelocities[corner8])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k3_y = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            xShift = 0.5 * spacing;
            yShift = 0.5 * spacing;
            // Do the same for stage 3 z velocities
            cellXCoordinate = floor((thirdXPosition - xShift) / spacing);
            cellYCoordinate = floor((thirdYPosition - yShift) / spacing);
            cellZCoordinate = floor(thirdZPosition / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (thirdXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = (thirdYPosition - yShift) - cellYCoordinate * spacing;
            deltaZ = thirdZPosition - cellZCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (1.0 - (deltaZ / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (1.0 - (deltaZ / spacing));
            w5 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w6 = (deltaX / spacing) * (1.0 - (deltaY / spacing)) * (deltaZ / spacing);
            w7 = (deltaX / spacing) * (deltaY / spacing) * (deltaZ / spacing);
            w8 = (1.0 - (deltaX / spacing)) * (deltaY / spacing) * (deltaZ / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
            corner2 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, cellZCoordinate);
            corner3 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner4 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), cellZCoordinate);
            corner5 = IX(cellXCoordinate, cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner6 = IX(fmin(cellXCoordinate + 1, gridLength - 2), cellYCoordinate, fmin(cellZCoordinate + 1, gridWidth - 2));
            corner7 = IX(fmin(cellXCoordinate + 1, gridLength - 2), fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            corner8 = IX(cellXCoordinate, fmin(cellYCoordinate + 1, gridHeight - 2), fmin(cellZCoordinate + 1, gridWidth - 2));
            // Proceed with transfer from the grid to particle
            offset = 1;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] == EMPTY || cellType[corner1 - offset] == EMPTY ? 0.0 : 1.0;
            valid2 = cellType[corner2] == EMPTY || cellType[corner2 - offset] == EMPTY ? 0.0 : 1.0;
            valid3 = cellType[corner3] == EMPTY || cellType[corner3 - offset] == EMPTY ? 0.0 : 1.0;
            valid4 = cellType[corner4] == EMPTY || cellType[corner4 - offset] == EMPTY ? 0.0 : 1.0;
            valid5 = cellType[corner5] == EMPTY || cellType[corner5 - offset] == EMPTY ? 0.0 : 1.0;
            valid6 = cellType[corner6] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid7 = cellType[corner7] == EMPTY || cellType[corner7 - offset] == EMPTY ? 0.0 : 1.0;
            valid8 = cellType[corner8] == EMPTY || cellType[corner8 - offset] == EMPTY ? 0.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4 + w5 * valid5 + w6 * valid6 + w7 * valid7 + w8 * valid8;
            velocity = particles[i].vz;
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * zVelocities[corner1] + valid2 * w2 * zVelocities[corner2] + valid3 * w3 * zVelocities[corner3] + valid4 * w4 * zVelocities[corner4] + valid5 * w5 * zVelocities[corner5] + valid6 * w6 * zVelocities[corner6] + valid7 * w7 * zVelocities[corner7] + valid8 * w8 * zVelocities[corner8]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (zVelocities[corner1] - prevZVelocities[corner1]) + valid2 * w2 * (zVelocities[corner2] - prevZVelocities[corner2]) + valid3 * w3 * (zVelocities[corner3] - prevZVelocities[corner3]) + valid4 * w4 * (zVelocities[corner4] - prevZVelocities[corner4]) + valid5 * w5 * (zVelocities[corner5] - prevZVelocities[corner5]) + valid6 * w6 * (zVelocities[corner6] - prevZVelocities[corner6]) + valid7 * w7 * (zVelocities[corner7] - prevZVelocities[corner7]) + valid8 * w8 * (zVelocities[corner8] - prevZVelocities[corner8])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k3_z = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            // Finally, update positions
            particles[i].x += ((2.0/9.0) * k1_x + (3.0/9.0) * k2_x + (5.0/9.0) * k3_x);
            particles[i].y += ((2.0/9.0) * k1_y + (3.0/9.0) * k2_y + (5.0/9.0) * k3_y);
            particles[i].z += ((2.0/9.0) * k1_z + (3.0/9.0) * k2_z + (5.0/9.0) * k3_z);

            // Also update velocities
            particles[i].vy += (-9.81 * dt);
        }
    }
    else    {
        for (int i = 0; i < numParticles; i++)  {
            // Only velocity in the Y direction is integrated - due to gravity
            particles[i].vy += (- dt * 9.81f);
            particles[i].x += (dt * particles[i].vx);
            particles[i].y += (dt * particles[i].vy);
            particles[i].z += (dt * particles[i].vz);
        }
    }
}

void Fluid::extrapolateVelocities() {
    // extend the velocity field throughout the grid (not just where the particles are)
    // use an empty vector to store the 'wavefront' for BFS
    std::vector<std::tuple<int, int, int>> wX;
    std::vector<std::tuple<int, int, int>> wY;
    std::vector<std::tuple<int, int, int>> wZ;
    // Part 1: X extrapolation
    // loop over the entire grid (barring solid cells at grid edges)
    for (int i = 1; i < gridLength; i++)    {
        for (int j = 1; j < gridHeight; j++) {
            for (int k = 1; k < gridWidth; k++) {
                if (cellType[IX(i, j, k)] == SOLID) continue;
                bool neighborXMarker = (xMarker[IX(i + 1, j, k)] == 0 || xMarker[IX(i - 1, j, k)] == 0 || xMarker[IX(i, j + 1, k)] == 0 || xMarker[IX(i, j - 1, k)] == 0 || xMarker[IX(i, j, k - 1)] == 0 || xMarker[IX(i, j, k + 1)] == 0);
                if (xMarker[IX(i, j, k)] != 0 && neighborXMarker)  {
                    xMarker[IX(i, j, k)] = 0;
                    wX.emplace_back(i, j, k);
                }
            }
        }
    }
    for (int t = 0; t < wX.size(); t++) {
        int i = std::get<0>(wX[t]);
        int j = std::get<1>(wX[t]);
        int k = std::get<2>(wX[t]);
        int numValidNeighbors = 0;
        int sum = 0;
        // check left neighbor
        if (cellType[IX(i - 1, j, k)] != SOLID) {
            if (xMarker[IX(i - 1, j, k)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i - 1, j, k)];
            }
            else if (xMarker[IX(i - 1, j, k)] == INT_MAX) {
                xMarker[IX(i - 1, j, k)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i - 1, j, k));
            }
        }
        // check right neighbor
        if (cellType[IX(i + 1, j, k)] != SOLID) {
            if (xMarker[IX(i + 1, j, k)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i + 1, j, k)];
            }
            else if (xMarker[IX(i + 1, j, k)] == INT_MAX)  {
                xMarker[IX(i + 1, j, k)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i + 1, j, k));
            }
        }
        // check bottom neighbor
        if (cellType[IX(i, j - 1, k)] != SOLID) {
            if (xMarker[IX(i, j - 1, k)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j - 1, k)];
            }
            else if (xMarker[IX(i, j - 1, k)] == INT_MAX)  {
                xMarker[IX(i, j - 1, k)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i, j - 1, k));
            }
        }
        // check top neighbor
        if (cellType[IX(i, j + 1, k)] != SOLID) {
            if (xMarker[IX(i, j + 1, k)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j + 1, k)];
            }
            else if (xMarker[IX(i, j + 1, k)] == INT_MAX)  {
                xMarker[IX(i, j + 1, k)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i, j + 1, k));
            }
        }
        // check back neighbor
        if (cellType[IX(i, j, k - 1)] != SOLID) {
            if (xMarker[IX(i, j, k - 1)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j, k - 1)];
            }
            else if (xMarker[IX(i, j, k - 1)] == INT_MAX)  {
                xMarker[IX(i, j, k - 1)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i, j, k - 1));
            }
        }
        // check front neighbor
        if (cellType[IX(i, j, k + 1)] != SOLID) {
            if (xMarker[IX(i, j, k + 1)] < xMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j, k + 1)];
            }
            else if (xMarker[IX(i, j, k + 1)] == INT_MAX)  {
                xMarker[IX(i, j, k + 1)] = xMarker[IX(i, j, k)] + 1;
                wX.push_back(std::tuple(i, j, k + 1));
            }
        }
        if (numValidNeighbors > 0) xVelocities[IX(i, j, k)] = sum / numValidNeighbors;
    }

    // Part 2: Y extrapolation
    // loop over the entire grid (barring solid cells at grid edges)
    for (int i = 1; i < gridLength; i++)    {
        for (int j = 1; j < gridHeight; j++) {
            for (int k = 1; k < gridWidth; k++) {
                if (cellType[IX(i, j, k)] == SOLID) continue;
                bool neighborYMarker = (yMarker[IX(i + 1, j, k)] == 0 || yMarker[IX(i - 1, j, k)] == 0 || yMarker[IX(i, j + 1, k)] == 0 || yMarker[IX(i, j - 1, k)] == 0 || yMarker[IX(i, j, k - 1) == 0 || yMarker[IX(i, j, k + 1)] == 0]);
                if (yMarker[IX(i, j, k)] != 0 && neighborYMarker)  {
                    yMarker[IX(i, j, k)] = 0;
                    wY.emplace_back(i, j, k);
                }
            }
        }
    }
    for (int t = 0; t < wY.size(); t++) {
        int i = std::get<0>(wY[t]);
        int j = std::get<1>(wY[t]);
        int k = std::get<2>(wY[t]);
        int numValidNeighbors = 0;
        int sum = 0;
        // check left neighbor
        if (cellType[IX(i - 1, j, k)] != SOLID) {
            if (yMarker[IX(i - 1, j, k)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i - 1, j, k)];
            } else if (yMarker[IX(i - 1, j, k)] == INT_MAX) {
                yMarker[IX(i - 1, j, k)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i - 1, j, k));
            }
        }
        // check right neighbor
        if (cellType[IX(i + 1, j, k)] != SOLID) {
            if (yMarker[IX(i + 1, j, k)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i + 1, j, k)];
            } else if (yMarker[IX(i + 1, j, k)] == INT_MAX) {
                yMarker[IX(i + 1, j, k)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i + 1, j, k));
            }
        }
        // check bottom neighbor
        if (cellType[IX(i, j - 1, k)] != SOLID) {
            if (yMarker[IX(i, j - 1, k)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j - 1, k)];
            } else if (yMarker[IX(i, j - 1, k)] == INT_MAX) {
                yMarker[IX(i, j - 1, k)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i, j - 1, k));
            }
        }
        // check top neighbor
        if (cellType[IX(i, j + 1, k)] != SOLID) {
            if (yMarker[IX(i, j + 1, k)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j + 1, k)];
            } else if (yMarker[IX(i, j + 1, k)] == INT_MAX) {
                yMarker[IX(i, j + 1, k)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i, j + 1, k));
            }
        }
        // check back neighbor
        if (cellType[IX(i, j, k - 1)] != SOLID) {
            if (yMarker[IX(i, j, k - 1)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j, k - 1)];
            } else if (yMarker[IX(i, j, k - 1)] == INT_MAX) {
                yMarker[IX(i, j, k - 1)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i, j, k - 1));
            }
        }
        // check front neighbor
        if (cellType[IX(i, j, k + 1)] != SOLID) {
            if (yMarker[IX(i, j, k + 1)] < yMarker[IX(i, j, k)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j, k + 1)];
            } else if (yMarker[IX(i, j, k + 1)] == INT_MAX) {
                yMarker[IX(i, j, k + 1)] = yMarker[IX(i, j, k)] + 1;
                wY.push_back(std::tuple(i, j, k + 1));
            }
        }
        if (numValidNeighbors > 0) yVelocities[IX(i, j, k)] = sum / numValidNeighbors;
    }

        // Part 3: Z extrapolation
        // loop over the entire grid (barring solid cells at grid edges)
        for (int i = 1; i < gridLength; i++) {
            for (int j = 1; j < gridHeight; j++) {
                for (int k = 1; k < gridWidth; k++) {
                    if (cellType[IX(i, j, k)] == SOLID) continue;
                    bool neighborZMarker = (zMarker[IX(i + 1, j, k)] == 0 || zMarker[IX(i - 1, j, k)] == 0 ||
                                            zMarker[IX(i, j + 1, k)] == 0 || zMarker[IX(i, j - 1, k)] == 0 ||
                                            zMarker[IX(i, j, k - 1) == 0 || zMarker[IX(i, j, k + 1)] == 0]);
                    if (zMarker[IX(i, j, k)] != 0 && neighborZMarker) {
                        zMarker[IX(i, j, k)] = 0;
                        wZ.emplace_back(i, j, k);
                    }
                }
            }
        }
        for (int t = 0; t < wZ.size(); t++) {
            int i = std::get<0>(wZ[t]);
            int j = std::get<1>(wZ[t]);
            int k = std::get<2>(wZ[t]);
            int numValidNeighbors = 0;
            int sum = 0;
            // check left neighbor
            if (cellType[IX(i - 1, j, k)] != SOLID) {
                if (zMarker[IX(i - 1, j, k)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i - 1, j, k)];
                } else if (zMarker[IX(i - 1, j, k)] == INT_MAX) {
                    zMarker[IX(i - 1, j, k)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i - 1, j, k));
                }
            }
            // check right neighbor
            if (cellType[IX(i + 1, j, k)] != SOLID) {
                if (zMarker[IX(i + 1, j, k)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i + 1, j, k)];
                } else if (zMarker[IX(i + 1, j, k)] == INT_MAX) {
                    zMarker[IX(i + 1, j, k)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i + 1, j, k));
                }
            }
            // check bottom neighbor
            if (cellType[IX(i, j - 1, k)] != SOLID) {
                if (zMarker[IX(i, j - 1, k)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i, j - 1, k)];
                } else if (zMarker[IX(i, j - 1, k)] == INT_MAX) {
                    zMarker[IX(i, j - 1, k)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i, j - 1, k));
                }
            }
            // check top neighbor
            if (cellType[IX(i, j + 1, k)] != SOLID) {
                if (zMarker[IX(i, j + 1, k)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i, j + 1, k)];
                } else if (zMarker[IX(i, j + 1, k)] == INT_MAX) {
                    zMarker[IX(i, j + 1, k)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i, j + 1, k));
                }
            }
            // check back neighbor
            if (cellType[IX(i, j, k - 1)] != SOLID) {
                if (zMarker[IX(i, j, k - 1)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i, j, k - 1)];
                } else if (zMarker[IX(i, j, k - 1)] == INT_MAX) {
                    zMarker[IX(i, j, k - 1)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i, j, k - 1));
                }
            }
            // check front neighbor
            if (cellType[IX(i, j, k + 1)] != SOLID) {
                if (zMarker[IX(i, j, k + 1)] < zMarker[IX(i, j, k)]) {
                    numValidNeighbors++;
                    sum += zVelocities[IX(i, j, k + 1)];
                } else if (zMarker[IX(i, j, k + 1)] == INT_MAX) {
                    zMarker[IX(i, j, k + 1)] = zMarker[IX(i, j, k)] + 1;
                    wZ.push_back(std::tuple(i, j, k + 1));
                }
            }
            if (numValidNeighbors > 0) zVelocities[IX(i, j, k)] = sum / numValidNeighbors;
        }
}