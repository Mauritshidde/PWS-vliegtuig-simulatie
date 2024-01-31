#include "MeshCube.h"

MeshCube::MeshCube(double vx, double vy, double vz, double setPressure, bool bound)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    boundary = bound;
    velocityX = 1;
    velocityY = 0.5;
    velocityZ = 0.1;
    pressureChanged = false;
}

MeshCube::~MeshCube()
{
}