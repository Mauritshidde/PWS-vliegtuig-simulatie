#include "MeshCube.h"

MeshCube::MeshCube(double vx, double vy, double vz, double setPressure, bool bound)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    boundary = bound;
    velocityX = 2;
    velocityY = 1;
    velocityZ = 1;
    pressureChanged = false;
}

MeshCube::~MeshCube()
{
}