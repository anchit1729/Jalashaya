//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

void Fluid::markFluidCells() {
    // Mark all boundary cells as solid cells to begin with
    for (int i = 0; i < gridLength; i++)  {
        for (int j = 0; j < gridHeight; j++)    {
            for (int k = 0; k < gridWidth; k++) {
                int cellIndex = IX(i, j, k);
                bool xBoundary = (i == 0 || i == gridLength - 1);
                bool yBoundary = (j == 0);//|| j == gridHeight - 1);
                bool zBoundary = (k == 0 || k == gridWidth - 1);
                if (xBoundary || yBoundary || zBoundary) cellType[cellIndex] = SOLID;
                else cellType[cellIndex] = EMPTY;
            }
        }
    }
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = Utils::clamp(floor(particles[i].x / spacing), 0, gridLength - 1);
        int cellYCoordinate = Utils::clamp(floor(particles[i].y / spacing), 0, gridHeight - 1);
        int cellZCoordinate = Utils::clamp(floor(particles[i].z / spacing), 0, gridWidth - 1);
        int cellIndex = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
        if (cellType[cellIndex] == EMPTY) cellType[cellIndex] = FLUID;
    }
}

void Fluid::transferVelocitiesToGrid() {
    // Transfer velocities of particles to their nearby grid side faces
    // First, as part of the FLIP method, we store a copy of all grid velocities
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    prevZVelocities = zVelocities;
    // Now, reset the velocity values
    std::fill(xVelocities.begin(), xVelocities.end(), 0.0);
    std::fill(yVelocities.begin(), yVelocities.end(), 0.0);
    std::fill(zVelocities.begin(), zVelocities.end(), 0.0);
    // Also reset the marker values
    std::fill(xMarker.begin(), xMarker.end(), INT_MAX);
    std::fill(yMarker.begin(), yMarker.end(), INT_MAX);
    std::fill(zMarker.begin(), zMarker.end(), INT_MAX);
    // Also, make sure delta values are set to 0.0
    std::fill(xR.begin(), xR.end(), 0.0);
    std::fill(yR.begin(), yR.end(), 0.0);
    std::fill(zR.begin(), zR.end(), 0.0);

    // Iterate over particles and mark fluid cells accordingly
//    for (int i = 0; i < numParticles; i++)  {
//        // Compute the cell coordinates
//        int cellXCoordinate = Utils::clamp(floor(particles[i].x / spacing), 0, gridLength - 1);
//        int cellYCoordinate = Utils::clamp(floor(particles[i].y / spacing), 0, gridHeight - 1);
//        int cellZCoordinate = Utils::clamp(floor(particles[i].z / spacing), 0, gridWidth - 1);
//        // Mark the cell as FLUID if not done already
//        int cellIndex = IX(cellXCoordinate, cellYCoordinate, cellZCoordinate);
//        if (cellType[cellIndex] == EMPTY) cellType[cellIndex] = FLUID;
//    }

    // Now, commence transfer to grid
    // Part 1: Transfer X velocities
    float yShift = 0.5 * spacing; // Shift downwards to account for staggered grid
    float zShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor(x / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor((y - yShift) / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor((z - zShift) / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = x - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        float deltaZ = (z - zShift) - cellZCoordinate * spacing;
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particles[i].vx;
        xVelocities[corner1] += particleVelocity * w1; xR[corner1] += w1; xMarker[corner1] = 0;
        xVelocities[corner2] += particleVelocity * w2; xR[corner2] += w2; xMarker[corner2] = 0;
        xVelocities[corner3] += particleVelocity * w3; xR[corner3] += w3; xMarker[corner3] = 0;
        xVelocities[corner4] += particleVelocity * w4; xR[corner4] += w4; xMarker[corner4] = 0;
        xVelocities[corner5] += particleVelocity * w5; xR[corner5] += w5; xMarker[corner5] = 0;
        xVelocities[corner6] += particleVelocity * w6; xR[corner6] += w6; xMarker[corner6] = 0;
        xVelocities[corner7] += particleVelocity * w7; xR[corner7] += w7; xMarker[corner7] = 0;
        xVelocities[corner8] += particleVelocity * w8; xR[corner8] += w8; xMarker[corner8] = 0;
    }
    // Wrap up the X velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < xVelocities.size(); i++)    {
        if (xR[i] > 0.0) xVelocities[i] /= xR[i];
    }

    // Part 2: Transfer Y velocities
    float xShift = 0.5 * spacing; // Shift sideways to account for staggered grid
    zShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor((x - xShift) / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor(y / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor((z - zShift) / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = y - cellYCoordinate * spacing;
        float deltaZ = (z - zShift) - cellZCoordinate * spacing;
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particles[i].vy;
        yVelocities[corner1] += particleVelocity * w1; yR[corner1] += w1; yMarker[corner1] = 0;
        yVelocities[corner2] += particleVelocity * w2; yR[corner2] += w2; yMarker[corner2] = 0;
        yVelocities[corner3] += particleVelocity * w3; yR[corner3] += w3; yMarker[corner3] = 0;
        yVelocities[corner4] += particleVelocity * w4; yR[corner4] += w4; yMarker[corner4] = 0;
        yVelocities[corner5] += particleVelocity * w5; yR[corner5] += w5; yMarker[corner5] = 0;
        yVelocities[corner6] += particleVelocity * w6; yR[corner6] += w6; yMarker[corner6] = 0;
        yVelocities[corner7] += particleVelocity * w7; yR[corner7] += w7; yMarker[corner7] = 0;
        yVelocities[corner8] += particleVelocity * w8; yR[corner8] += w8; yMarker[corner8] = 0;
    }
    // Wrap up the Y velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < yVelocities.size(); i++)    {
        if (yR[i] > 0.0) yVelocities[i] /= yR[i];
    }

    // Part 3: Transfer Z velocities
    xShift = 0.5 * spacing; // Shift sideways to account for staggered grid
    yShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
        float x = Utils::clamp(particles[i].x, spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particles[i].y, spacing, (gridHeight - 1) * spacing);
        float z = Utils::clamp(particles[i].z, spacing, (gridWidth - 1) * spacing);
        int cellXCoordinate = fmin(floor((x - xShift) / spacing), gridLength - 2);
        int cellYCoordinate = fmin(floor((y - yShift) / spacing), gridHeight - 2);
        int cellZCoordinate = fmin(floor(z / spacing), gridWidth - 2);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        float deltaZ = z - cellZCoordinate * spacing;
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particles[i].vz;
        zVelocities[corner1] += particleVelocity * w1; zR[corner1] += w1; zMarker[corner1] = 0;
        zVelocities[corner2] += particleVelocity * w2; zR[corner2] += w2; zMarker[corner2] = 0;
        zVelocities[corner3] += particleVelocity * w3; zR[corner3] += w3; zMarker[corner3] = 0;
        zVelocities[corner4] += particleVelocity * w4; zR[corner4] += w4; zMarker[corner4] = 0;
        zVelocities[corner5] += particleVelocity * w5; zR[corner5] += w5; zMarker[corner5] = 0;
        zVelocities[corner6] += particleVelocity * w6; zR[corner6] += w6; zMarker[corner6] = 0;
        zVelocities[corner7] += particleVelocity * w7; zR[corner7] += w7; zMarker[corner7] = 0;
        zVelocities[corner8] += particleVelocity * w8; zR[corner8] += w8; zMarker[corner8] = 0;
    }
    // Wrap up the Z velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < zVelocities.size(); i++)    {
        if (zR[i] > 0.0) zVelocities[i] /= zR[i];
    }

    // restore grid velocities for solid cells
    for (int i = 0; i < gridLength; i++) {
        for (int j = 0; j < gridHeight; j++) {
            for (int k = 0; k < gridWidth; k++) {
                if (cellType[IX(i, j, k)] == SOLID || (i > 0 && cellType[IX(i - 1, j, k)] == SOLID))   {
                    xVelocities[IX(i, j, k)] = prevXVelocities[IX(i, j, k)];
                }
                if (cellType[IX(i, j, k)] == SOLID || (j > 0 && cellType[IX(i, j - 1, k)] == SOLID))   {
                    yVelocities[IX(i, j, k)] = prevYVelocities[IX(i, j, k)];
                }
                if (cellType[IX(i, j, k)] == SOLID || (k > 0 && cellType[IX(i, j, k - 1)] == SOLID))   {
                    zVelocities[IX(i, j, k)] = prevZVelocities[IX(i, j, k)];
                }
            }
        }
    }
}