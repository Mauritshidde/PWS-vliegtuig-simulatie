#include "MeshCube.h"

MeshCube::MeshCube(float vx, float vy, float vz, float setPressure, bool bound)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    boundary = bound;
}

MeshCube::~MeshCube()
{
}