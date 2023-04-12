#include <iostream>
#include "Fluid.h"
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include <fstream>

int main() {
    std::cout << "2D PIC/FLIP Simulation.\n" << std::endl;
    Fluid fluid;
    polyscope::init();
    polyscope::PointCloud* psCloud;
    std::vector<glm::vec3> points;

    for (int i = 0; i < fluid.numParticles; i++)    {
        float x = fluid.particles[i].x;
        float y = fluid.particles[i].y;
        float z = fluid.particles[i].z;
        points.push_back(glm::vec3{x, y, z});
    }


    psCloud = polyscope::registerPointCloud("Flip Fluid", points);
    psCloud->setPointRadius(fluid.particleRadius*0.001);
    //psCloud->setPointRenderMode(polyscope::PointRenderMode::Quad);
    // simulate 480 frames total
    int frameNumber = 0;
    polyscope::show();
    for (int iter = 0; iter < 1200; iter++)    {
        std::ofstream csv_file("coordinates/frame" + std::to_string(iter+1) + ".txt");
        for (int i = 0; i < fluid.numParticles; i++)    {
            float x = fluid.particles[i].x;
            float y = fluid.particles[i].y;
            float z = fluid.particles[i].z;
            float vx = fluid.particles[i].vx;
            float vy = fluid.particles[i].vy;
            float vz = fluid.particles[i].vz;
            points[i] = glm::vec3{x, y, z};
            csv_file << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << "\n";
        }
        if (frameNumber == -1)  {
            polyscope::show();
        }
        std::string frameName = "./screenshots/Frame" + std::to_string(frameNumber) + ".png";
        psCloud->updatePointPositions(points);
        polyscope::screenshot(frameName);
        frameNumber++;
        std::cout << "Rendered " << frameNumber << " frames\n";
        fluid.simulateFluid();
        //polyscope::show();
        csv_file.close();
    }
    polyscope::shutdown();


    return 0;
}
