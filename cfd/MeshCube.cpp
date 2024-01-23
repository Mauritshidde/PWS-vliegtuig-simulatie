#include "MeshCube.h"

MeshCube::MeshCube(float vx, float vy, float vz, float setPressure, float setDensity, bool bound)
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