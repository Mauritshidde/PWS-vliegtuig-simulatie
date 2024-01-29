#include "MeshCube.h"

MeshCube::MeshCube(double vx, double vy, double vz, double setPressure, bool bound)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    boundary = bound;
    velocityX = 1;
    velocityY = 0;
    velocityZ = 0;
    pressureChanged = false;
}

MeshCube::~MeshCube()
{
}