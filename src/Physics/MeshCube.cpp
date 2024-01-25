#include "MeshCube.h"

MeshCube::MeshCube(double vx, double vy, double vz, double setPressure, double setDensity, bool bound)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    density = setDensity;
    boundary = bound;
}

MeshCube::~MeshCube()
{
}