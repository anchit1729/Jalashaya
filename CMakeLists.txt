cmake_minimum_required(VERSION 3.24)
project(FluidSim_3D)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(polyscope)
add_executable(FluidSim_3D main.cpp Fluid.h Advection.cpp Projection.cpp SpatialHashing.cpp TransferToGrid.cpp TransferToParticles.cpp Utils.h Particle.h Fluid.cpp)
target_link_libraries(FluidSim_3D polyscope)